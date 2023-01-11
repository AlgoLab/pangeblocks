"""Solve the ILP formulation:
Given an MSA, and a set of blocks such that the union of them covers all the MSA,
Find a non-overlapping set the blocks covering the entire MSA
# """
from collections import deque
from itertools import chain
from sys import getsizeof, stderr
import time
from collections import defaultdict
from tqdm import tqdm
from ..blocks import Block
from pathlib import Path
from Bio import AlignIO
from gurobipy import GRB
import gurobipy as gp
from src.blocks.analyzer import BlockAnalyzer
import logging
logging.basicConfig(level=logging.ERROR)
# log = Logger(name="opt", level="DEBUG")

try:
    from reprlib import repr
except ImportError:
    pass


def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}

    """
    def dict_handler(d): return chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                    }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    # estimate sizeof object without __sizeof__
    default_size = getsizeof(0)

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        if verbose:
            print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)



class Optimization:
    def __init__(self, blocks, path_msa, path_save_ilp=None, log_level=logging.ERROR):

        self.input_blocks = blocks
        # Each block is a tuple (K, i, j, label) where
        # K is a tuple of rows,
        # i and j are the first and last column
        # label is the string of the block
        msa, n_seqs, n_cols = self.load_msa(path_msa)
        self.msa = msa
        self.n_seqs = n_seqs
        self.n_cols = n_cols
        self.path_save_ilp = path_save_ilp
        logging.getLogger().setLevel(log_level)

    def __call__(self, return_times: bool = False):
        "Solve ILP formulation"
        times = dict()
        ti = time.time()

        # covering_by_position is a dictionary with key (r,c) and value the list
        # of indices of the blocks that include the position (r,c)
        #
        # to speed up the process, we iterate over the blocks and append the
        # block to all the positions it covers
        
        # save vertical blocks and their covered positions in the MSA
        vertical_blocks = []
        covered_positions = defaultdict(list)

        # save positions not covered by vertical blocks separatedly
        msa_positions = []
        covering_by_position = defaultdict(list)
        for idx, block in enumerate(self.input_blocks):
            
            # vertical blocks: will be chosen in the optimal solution
            if len(block.K) == self.n_seqs: 
                vertical_blocks.append(idx)
            
            # all other blocks
            for r in block.K:
                for c in range(block.i, block.j + 1):
                    covering_by_position[(r, c)].append(idx)

        # write idx for MSA positions (row, col)
        msa_positions = [(r, c) for r in range(self.n_seqs)
                         for c in range(self.n_cols)]


        tf = time.time()
        times["init"] = round(tf - ti, 3)

        # Create the model
        ti = time.time()
        model = gp.Model("pangeblocks")

        # Threads
        model.setParam(GRB.Param.Threads, 8)

        # define variables
        # C(b) = 1 if block b is selected
        # U(r,c) = 1 if position (r,c) is covered by at least one block
        C = model.addVars(range(len(self.input_blocks)),
                          vtype=GRB.BINARY, name="C")
        for block in range(len(self.input_blocks)):
            logging.info(
                f"variable:C({block}) = {self.input_blocks[block].str()}")
        for pos in msa_positions:
            logging.info(f"variable:U({pos})")
        U = model.addVars(msa_positions, vtype=GRB.BINARY, name="U")

        # Constraints
        for r, c in tqdm(msa_positions):

            blocks_rc = covering_by_position[(r, c)]
            
            if len(blocks_rc) > 0:
                # 1. U[r,c] = 1 implies that at least one block covers the position
                model.addConstr(
                    U[r, c] <= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint1({r},{c})"
                )

                # 2. each position in the MSA is covered at most by one block
                model.addConstr(
                    1 >= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint2({r},{c})"
                )
                logging.debug(f"constraint2({r},{c}) covered by {blocks_rc}")

                # 3. each position of the MSA is covered AT LEAST by one block
                model.addConstr(U[r, c] >= 1, name=f"constraint3({r},{c})")
                logging.debug(f"constraint3({r},{c})")

        tf = time.time()
        times["constraints1-2-3"] = round(tf - ti, 3)
        
        # constraint 4: vertical blocks are part of the solution
        for idx in vertical_blocks:
            model.addConstr(C[idx] == 1)
            logging.info(
                f"constraint4: vertical block ({idx}) - {self.input_blocks[idx].str()}"
            )

        ti = time.time()
        # TODO: include input to decide which objective function to use
        # Objective function
        # minimize the number of blocks (nodes)
        # model.setObjective(C.sum("*", "*", "*"), GRB.MINIMIZE)
        
        # minimize the total length of the graph (number of characters)
        # model.setObjective(
        #     sum(len(block.label)*C[idx] for idx, block in enumerate(self.input_blocks)), 
        #     GRB.MINIMIZE
        # )

        # minimize the number of blocks penalizing shorter blocks
        MIN_LEN = 2 # penalize blocks with label less than MIN_LEN
        PENALIZATION = 3 # costly than the other ones
        model.setObjective(
            sum( 
                (PENALIZATION if len(block.label)<MIN_LEN else 1)*C[idx] 
                for idx, block in enumerate(self.input_blocks)
                ), 
            GRB.MINIMIZE
        )

        logging.info("Begin ILP")

        model.optimize()
        logging.info("End ILP")
        tf = time.time()
        times["optimization"] = round(tf - ti, 3)

        if self.path_save_ilp:
            Path(self.path_save_ilp).parent.mkdir(exist_ok=True, parents=True)
            model.write(self.path_save_ilp)

        try:
            solution_C = model.getAttr("X", C)
        except:
            raise ("No solution")

        # filter optimal coverage of blocks for the MSA
        ti = time.time()
        optimal_coverage = []
        for k, v in solution_C.items():

            if v > 0:
                logging.info(
                    f"Optimal Solution: {k}, {self.input_blocks[k].str()}")
                optimal_coverage.append(self.input_blocks[k])
        tf = time.time()
        times["solution as blocks"] = round(tf - ti, 3)

        if return_times is True:
            return optimal_coverage, times
        return optimal_coverage

    def load_msa(self, path_msa):
        "return alignment, number of sequences and columns"
        # load MSA
        align = AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)

        return align, n_seqs, n_cols