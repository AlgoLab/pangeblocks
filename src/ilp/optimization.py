"""Solve the ILP formulation:
Given an MSA, and a set of blocks such that the union of them covers all the MSA,
Find a non-overlapping set the blocks covering the entire MSA
# """
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
        covering_by_position = defaultdict(list)
        for idx, block in enumerate(self.input_blocks):
            for r in block.K:
                logging.debug(f"block: {block.str()}/{set(block.K)}")
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
        U = model.addVars(msa_positions, vtype=GRB.BINARY, name="U")
        for pos in msa_positions:
            logging.info(f"variable:U({pos})")

        # Constraints
        for r, c in tqdm(msa_positions):
            blocks_rc = covering_by_position[(r, c)]
            if len(blocks_rc) > 0:
                # 1. each position in the MSA is covered ONLY ONCE
                model.addConstr(
                    U[r, c] <= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint1({r},{c})"
                )
                logging.info(f"constraint1({r},{c}) covered by {blocks_rc}")

                # 2. each position of the MSA is covered AT LEAST by one block
                model.addConstr(U[r, c] >= 1, name=f"constraint2({r},{c})")
                logging.info(f"constraint2({r},{c})")

        tf = time.time()
        times["constraints1-2"] = round(tf - ti, 3)

        # 3. overlapping blocks cannot be chosen
        # sort all blocks,
        ti = time.time()

        intersection = self._list_inter_blocks(self.input_blocks)

        for idx1, idx2 in tqdm(intersection):
            # if the blocks intersect, then create the restriction
            block1 = self.input_blocks[idx1]
            block2 = self.input_blocks[idx2]
            name_constraint = f"constraint3({idx1})-({idx2})"
            model.addConstr(C[idx1] + C[idx2] <= 1, name=name_constraint)
            logging.info(
                f"constraint3({idx1})-({idx2}), {block1.str()} - {block2.str()}")

        tf = time.time()
        times["constraint3"] = round(tf - ti, 3)

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

            # id_block, i, j = k
            if v > 0:
                logging.info(
                    f"Optimal Solution: {k}, {self.input_blocks[k].str()}")
                optimal_coverage.append(self.input_blocks[k])
        #         K = (int(seq) for seq in id_block_to_K[id_block].split(","))
        #         label = id_block_to_labels[id_block]
        #         optimal_coverage.append(Block(K, i, j, label))
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

    def get_common_cols(self, block1, block2):
        if block1.j < block2.i:
            return []
        intervals = [(block1.i, block1.j), (block2.i, block2.j)]
        start, end = intervals.pop()
        while intervals:
            start_temp, end_temp = intervals.pop()
            start = max(start, start_temp)
            end = min(end, end_temp)
        return [start, end]

    def _list_inter_blocks(self, list_blocks: list[Block]) -> list[tuple]:
        """returns a list of tuples(idx_block1, idx_block2) with all pairs
        of blocks that intersects"""
        # "list of indexes (in a sorted list by i) of pairs of blocks with non-empty intersection"

        # sort blocks by starting position and save the order in the original list
        sorted_blocks = sorted(enumerate(self.input_blocks),
                               key=lambda block: block[1].i)

        # save pairs of indexes for the sorted blocks that intersect
        # (pos1,pos2): positions in the sorted list
        # (orig_pos1, orig_pos2): positions in the order of the input list
        intersections = []
        visited_blocks = defaultdict(set)
        current_block_idx = 0
        blocks_ending = defaultdict(set)

        # For each column, we will compute the index of blocks that
        # start in that column, and store those indices in blocks_by_start
        block_by_start = defaultdict(list)
        for idx, block in sorted_blocks:
            block_by_start[block.i].append(idx)

        # We will iterate over the columns of the MSA.
        # For each column, we will
        # (1) compute the intersections between blocks starting in that column
        # and the blocks in visited_blocks
        # (2) remove from visited_blocks the blocks that end in the current column
        for column in range(self.n_cols):
            logging.debug(f"current column: {column}")
            logging.debug(f"visited blocks: {visited_blocks}")
            # notice that in block_by_start[column] we have the blocks that
            # start in the current column.
            # We are going to save them in the current_blocks set
            current_blocks = dict([(idx, self.input_blocks[idx])
                                  for idx in block_by_start[column]])
            logging.debug(f"current blocks: {current_blocks}")

            # (1) compute the intersections between blocks in block_by_start[column] and
            # blocks in visited_blocks
            # Notice that all visited blocks span the current column, so
            # there is no need to check if they intersect in the current column
            for idx1, block1 in current_blocks.items():
                k1 = set(block1.K)
                logging.debug(f"current block: {idx1}, {block1.str()}")
                for idx2, block2 in visited_blocks.items():
                    logging.debug(f"first block: {idx1}, {block1.str()}")
                    logging.debug(f"second block: {idx2}, {block2.str()}")
                    logging.debug(f"sets: {k1}, {set(block2.K)}")
                    if not k1.isdisjoint(set(block2.K)):
                        intersections.append((idx1, idx2))
                        logging.info(
                            f"intersect: {idx1}, {block1.str()}, {idx2}, {block2.str()}")

            # Add all blocks in current_blocks to the set of visited blocks
            visited_blocks = visited_blocks | current_blocks
            for idx1, current_block in current_blocks.items():
                blocks_ending[current_block.j].add(idx1)

            #  Now that we have processed all the blocks that start in the
            #  current column, we can remove the blocks that end in the current
            #  column from the list of visited blocks
            for idx in blocks_ending[column]:
                del visited_blocks[idx]

        return intersections
