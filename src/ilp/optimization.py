"""Solve the ILP formulation:
Given an MSA, and a set of blocks such that the union of them covers all the MSA,
Find a non-overlapping set the blocks covering the entire MSA
# """
from collections import deque
from itertools import chain
from sys import getsizeof, stderr
import time
from collections import defaultdict
from ..blocks import Block
from pathlib import Path
from Bio import AlignIO
from gurobipy import GRB
import gurobipy as gp
from src.blocks.analyzer import BlockAnalyzer
import logging
logging.basicConfig(level=logging.ERROR,
                    format='%(asctime)s. %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')
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

        # all blocks we receive involve at least two sequences, but we need to
        # add dummy blocks for the cell of the MSA.
        # This allows to guarantee that the ILP always has a feasible solution
        dummy_blocks = [Block((r,), c, c, self.msa[r,c]) for r in range(self.n_seqs)
                        for c in range(self.n_cols)]
        all_blocks = self.input_blocks + dummy_blocks
        # covering_by_position is a dictionary with key (r,c) and value the list
        # of indices of the blocks that include the position (r,c)
        #
        # to speed up the process, we iterate over the blocks and append the
        # block to all the positions it covers
        covering_by_position = defaultdict(list)
        left_maximal_vertical_blocks = {}

        # We keep only maximal vertical blocks and we achieve that in two
        # phases:
        # 1. for each beginning column, we keep the block with the largest
        #    end column
        # 2. we do a second scan of vertical blocks and, for each end
        #    column, we keep the block with the smallest beginning column
        for idx, block in enumerate(all_blocks):
            logging.debug(f"input block: {block.str()}")
            if len(block.K) == self.n_seqs:
                if not block.i in left_maximal_vertical_blocks or block.j > left_maximal_vertical_blocks[block.i]["end"]:
                    left_maximal_vertical_blocks[block.i] = {
                        "idx": idx, "end": block.j}
        vertical_blocks = {}
        for begin, block in left_maximal_vertical_blocks.items():
            if not block["end"] in vertical_blocks or begin < vertical_blocks[block["end"]]["begin"]:
                vertical_blocks[block["end"]] = {
                    "idx": block["idx"], "begin": begin}
        #
        # We keep a set of all columns that are covered by a vertical block:
        # those columns will not be involved in the U[r,c] variables and in the
        # corresponding covering constraints
        covered_by_vertical_block = set()
        for begin, item in vertical_blocks.items():
            block = all_blocks[item['idx']]
            logging.info(
                f"Vertical block: {block.str()}")
            for col in range(block.i, block.j + 1):
                covered_by_vertical_block.add(col)
        logging.info(
            f"Covered by vertical blocks: {covered_by_vertical_block}")
        logging.info(
            f"No. covered by vertical blocks: {len(covered_by_vertical_block)} out of {self.n_cols}")

        # We compute the set disjoint_vertical, corresponding to the
        # blocks that are disjoint from vertical blocks.
        # The blocks that we will encode with a C variable in the
        # ILP correspond to the blocks that are disjoint from vertical blocks or
        # are vertical blocks themselves.
        disjoint_vertical = set()
        for idx, block in enumerate(all_blocks):
            if block.i in covered_by_vertical_block or block.j in covered_by_vertical_block:
                pass
            if set(range(block.i, block.j + 1)).isdisjoint(covered_by_vertical_block):
                disjoint_vertical.add(idx)
                logging.debug(f"disjoint vertical block: {block.str()}")

        c_variables = list(disjoint_vertical) + [item["idx"]
                                                 for item in vertical_blocks.values()]
        # msa_positions is a list of all positions (r,c) that are required to be
        # covered. We exclude the positions covered by vertical blocks, since
        # they will be guaranteed to be covered, as an effect of the fact that
        # the corresponding C variables will be set to 1
        msa_positions = [(r, c) for r in range(self.n_seqs)
                         for c in set(range(self.n_cols)) - covered_by_vertical_block]

        # covering_by_position is a dictionary with key (r,c) and value the list
        # of indices of the blocks that cover the position (r,c)
        for idx in c_variables:
            block = all_blocks[idx]
            logging.debug(f"block: {block.str()}")
            for r in block.K:
                for c in range(block.i, block.j + 1):
                    covering_by_position[(r, c)].append(idx)

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
        # S(r,c) = 1 if position (r,c) is covered by a 1-cell block
        C = model.addVars(c_variables,
                          vtype=GRB.BINARY, name="C")
        for block in c_variables:
            logging.info(
                f"variable:C({block}) = {all_blocks[block].str()}")
        U = model.addVars(msa_positions, vtype=GRB.BINARY, name="U")
        for pos in msa_positions:
            logging.info(f"variable:U({pos})")

        # Constraints
        # All C(b) variables corresponding to vertical blocks are set to 1
        for end, item in vertical_blocks.items():
            model.addConstr(C[item["idx"]] == 1,
                            name=f"vertical_constraint({item['idx']})")

        for r, c in msa_positions:
            blocks_rc = covering_by_position[(r, c)]
            if len(blocks_rc) > 0:
                # U[r,c] = 1 implies that at least one block covers the position
                model.addConstr(
                    U[r, c] <= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint1({r},{c})"
                )
                # 1. each position in the MSA is covered by at most one block
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
        ti = time.time()

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
                    f"Optimal Solution: {k}, {all_blocks[k].str()}")
                optimal_coverage.append(all_blocks[k])
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
