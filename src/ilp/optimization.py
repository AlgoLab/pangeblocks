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
                # # U[r,c] = 1 implies that at least one block covers the position
                model.addConstr(
                    U[r, c] <= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint1({r},{c})"
                )

                # 1. each position in the MSA is covered at most by one block
                model.addConstr(
                    1 >= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint1({r},{c})"
                )
                logging.debug(f"constraint1({r},{c}) covered by {blocks_rc}")

                # 2. each position of the MSA is covered AT LEAST by one block
                model.addConstr(U[r, c] >= 1, name=f"constraint2({r},{c})")
                logging.debug(f"constraint2({r},{c})")

        tf = time.time()
        times["constraints1-2"] = round(tf - ti, 3)

        # 3. overlapping blocks cannot be chosen
        # sort all blocks,
        # ti = time.time()

        # logging.info("Constraints 2")
        # for idx1, idx2 in (self._list_inter_blocks(self.input_blocks)):
        #     # if the blocks intersect, then create the restriction
        #     block1 = self.input_blocks[idx1]
        #     block2 = self.input_blocks[idx2]
        #     name_constraint = f"constraint3({idx1})-({idx2})"
        #     model.addConstr(C[idx1] + C[idx2] <= 1, name=name_constraint)
        #     logging.debug(
        #         f"constraint3({idx1})-({idx2}), {block1.str()} - {block2.str()}")

        # tf = time.time()
        # times["constraint3"] = round(tf - ti, 3)
        
        # constraint 4: vertical blocks are part of the solution
        for idx in vertical_blocks:
            model.addConstr(C[idx] == 1)
            logging.info(
                f"constraint4: vertical block ({idx}) - {self.input_blocks[idx].str()}"
            )

        ti = time.time()
        # Objective function
        model.setObjective(C.sum("*", "*", "*"), GRB.MINIMIZE)
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

    # def _list_inter_blocks(self, list_blocks: list[Block]) -> list[tuple]:
    #     """returns a the tuples(idx_block1, idx_block2) for each pair
    #     of blocks that intersects.
    #     It is implemented as a generator, since the list might grow too much."""


    #     # We organize the blocks in two dictionaries:
    #     # blocks_ending: for each column, the blocks that end in that column
    #     # block_by_start: for each column, the blocks that start in that column
    #     blocks_ending = defaultdict(set)
    #     block_by_start = defaultdict(list)
    #     block_by_end = defaultdict(list)
    #     for idx, block in self.input_blocks:
    #         block_by_start[block.i].append(idx)
    #         block_by_end[block.j].append(idx)
    #     for column in block_by_end.keys():
    #         blocks_ending[column] = set(block_by_end[column])

    #     # We will iterate over the columns of the MSA.
    #     # We maintain a set of blocks that span the current column
    #     # (current_blocks) and a set of blocks that have been visited previously
    #     # (visited_blocks) and that span the current column.
    #     # For each column, we will
    #     # (1) compute the intersections between blocks starting in that column
    #     # and the blocks in visited_blocks
    #     # (2) remove from visited_blocks the blocks that end in the current
    #     # column (so that the next iteration will not consider them)
    #     visited_blocks = set()
    #     for column in range(self.n_cols):
    #         logging.info(f"current column: {column}")
    #         # notice that in block_by_start[column] we have the blocks that
    #         # start in the current column.
    #         # We are going to save them in the current_blocks set
    #         current_blocks = set([(idx, self.input_blocks[idx])
    #                               for idx in block_by_start[column]])
    #         logging.debug(f"current blocks: {current_blocks}")

    #         # (1) compute the intersections between blocks in block_by_start[column] and
    #         # blocks in visited_blocks
    #         # Notice that all visited blocks span the current column, so
    #         # there is no need to check if they intersect in the current column
    #         for idx1, block1 in current_blocks:
    #             k1 = set(block1.K)
    #             logging.debug(f"current block: {idx1}, {block1.str()}")
    #             for idx2, block2 in visited_blocks:
    #                 logging.debug(f"first block: {idx1}, {block1.str()}")
    #                 logging.debug(f"second block: {idx2}, {block2.str()}")
    #                 logging.debug(f"sets: {k1}, {set(block2.K)}")
    #                 if not k1.isdisjoint(set(block2.K)):
    #                     yield (idx1, idx2)
    #                     logging.debug(
    #                         f"intersect: {idx1}, {block1.str()}, {idx2}, {block2.str()}")
    #         logging.debug(
    #             f"Size visited_blocks: {total_size(visited_blocks):_}, len: {len(visited_blocks)}")
    #         logging.debug(
    #             f"Size blocks_ending: {total_size(blocks_ending):_}, len: {len(blocks_ending)}")
    #         logging.debug(
    #             f"Size block_by_start: {total_size(block_by_start):_}, len: {len(block_by_start)}")
    #         # Add all blocks in current_blocks to the set of visited blocks
    #         visited_blocks |= current_blocks
    #         #  Now that we have processed all the blocks that start in the
    #         #  current column, we can remove the blocks that end in the current
    #         #  column from the list of visited blocks
    #         visited_blocks -= blocks_ending[column]