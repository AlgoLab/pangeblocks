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

class Optimization:
    def __init__(self, blocks, path_msa, path_save_ilp=None, log_level=logging.ERROR, **kwargs):

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
        self.obj_function = kwargs.get("obj_function","nodes")
        self.penalization = kwargs.get("penalization",1)
        self.min_len = kwargs.get("min_len",1)
        self.time_limit = kwargs.get("time_limit", 180)
        logging.getLogger().setLevel(log_level)

    def __call__(self, return_times: bool = False):
        "Solve ILP formulation"
        times = dict()
        ti = time.time()

        all_blocks = self.input_blocks
        n_blocks = len(all_blocks)
        
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
        private_blocks = []
        for idx, block in enumerate(all_blocks):
            logging.info(f"analyzing {idx} out of {n_blocks}")
            if set(range(block.i, block.j + 1)).isdisjoint(covered_by_vertical_block):
                # The current block is disjoint from vertical blocks
                disjoint_vertical.add(idx)
                logging.debug(f"disjoint vertical block: {block.str()}")
            elif not set(range(block.i, block.j + 1)).issubset(covered_by_vertical_block):
                # The current block overlaps the vertical blocks.
                # We need to check if the block intersects with at least two
                # vertical blocks: in that case we need to manually decompose
                # it.

                # private_cols is the set of columns of the current block that
                # are not in any vertical blocks.
                # We need to extract from those columns the regions with
                # consecutive columns, so that we can decompose the current
                # block into the sub-blocks divided by the vertical blocks.
                private_cols = set(range(block.i, block.j + 1)) - covered_by_vertical_block
                sorted_cols = sorted(private_cols)
                private_regions = []
                begin = -1
                for idx, col in enumerate(sorted_cols):
                    if begin < 0:
                        begin = col
                    if idx == len(sorted_cols) - 1 or col < sorted_cols[idx + 1] - 1:
                        private_regions.append((begin, col))
                        begin = -1

                num_intersecting_vertical = len(private_regions) - 1
                if sorted_cols[0] in covered_by_vertical_block:
                    num_intersecting_vertical += 1
                if sorted_cols[-1] in covered_by_vertical_block:
                    num_intersecting_vertical += 1
                if num_intersecting_vertical >= 2:
                    logging.debug(
                        f"block {block.str()} intersects with at least two vertical blocks")
                    logging.debug(f"sorted_cols: {sorted_cols}")
                    logging.debug(f"private_regions: {private_regions}")
                    for begin, end in private_regions:
                        label = "".join([self.msa[block.K[0], i] for i in range(begin, end + 1)])
                        new_block = Block(block.K, begin, end, label)
                        private_blocks.append(new_block)
                        logging.debug(f"Adding private block: {new_block.str()} to {private_blocks}")
            
        first_private_block = len(all_blocks)
        # all_blocks += enumerate(private_blocks, first_private_block)
        all_blocks += private_blocks

        c_variables = list(disjoint_vertical) + [item["idx"]
                                                 for item in vertical_blocks.values()] + list(range(first_private_block, len(all_blocks)))
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
        model.setParam(GRB.Param.Threads, 16)

        # Time Limit
        model.setParam(GRB.Param.TimeLimit, float(60*self.time_limit))
        
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

        # # constraint 4: vertical blocks are part of the solution
        # for idx in vertical_blocks:
        #     model.addConstr(C[idx] == 1)
        #     logging.info(
        #         f"constraint4: vertical block ({idx}) - {self.input_blocks[idx].str()}"
        #     )

        ti = time.time()
        # TODO: include input to decide which objective function to use
        # Objective function
        if self.obj_function == "nodes":
            # minimize the number of blocks (nodes)
            model.setObjective(C.sum("*", "*", "*"), GRB.MINIMIZE)
        elif self.obj_function == "strings":
            # minimize the total length of the graph (number of characters)
            model.setObjective(
                sum(
                    len(all_blocks[idx].label)*C[idx] for idx in c_variables
                    ), 
                GRB.MINIMIZE
            )
        elif self.obj_function == "weighted":
            # minimize the number of blocks penalizing shorter blocks
            MIN_LEN = self.min_len # penalize blocks with label less than MIN_LEN
            PENALIZATION = self.penalization # costly than the other ones
            model.setObjective(
                sum( 
                    (PENALIZATION if len(all_blocks[idx].label)<MIN_LEN else 1)*C[idx] 
                    for idx in c_variables
                    ),
                GRB.MINIMIZE
            )

        logging.info(f"Begin ILP with Objective function: {self.obj_function}")
        if self.obj_function == "weighted":
            logging.info(f"penalization: {self.penalization}")
            logging.info(f"minimum length: {self.min_len}")

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