"""Solve the ILP formulation:
Given an MSA, and a set of blocks such that the union of them covers all the MSA,
Find a non-overlapping set the blocks covering the entire MSA
# """
from collections import deque
from itertools import chain
from sys import getsizeof, stderr
import time
from collections import defaultdict
from blocks import Block
from pathlib import Path
from Bio import AlignIO
from gurobipy import GRB, LinExpr
import gurobipy as gp
from blocks import BlockAnalyzer
from .losses import (
    loss_nodes,
    loss_strings,
    loss_weighted,
    loss_depth
)
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
        self.obj_function = kwargs.get("obj_function", "nodes")
        self.penalization = kwargs.get("penalization", 1)
        self.min_len = kwargs.get("min_len", 1)
        self.min_coverage = kwargs.get("min_coverage", 1)
        self.time_limit = kwargs.get("time_limit", 180)
        logging.getLogger().setLevel(log_level)

    def __call__(self, return_times: bool = False):
        "Solve ILP formulation"
        times = dict()
        ti = time.time()

        n_blocks = len(self.input_blocks)
        logging.info("number of blocks %s", n_blocks)
        for idx, block in enumerate(self.input_blocks):
            logging.debug("input block %s: %s" % (idx, block.str()))

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
        logging.info("loop: left maximal vertical blocks")
        for idx, block in enumerate(self.input_blocks):
            # logging.debug("input block: %s", block.str())
            if len(block.K) == self.n_seqs:
                if not block.start in left_maximal_vertical_blocks or block.end > left_maximal_vertical_blocks[block.start]["end"]:
                    left_maximal_vertical_blocks[block.start] = {
                        "idx": idx, "end": block.end}

        logging.info("loop: vertical blocks")
        vertical_blocks = {}
        for begin, block in left_maximal_vertical_blocks.items():
            if not block["end"] in vertical_blocks or begin < vertical_blocks[block["end"]]["begin"]:
                vertical_blocks[block["end"]] = {
                    "idx": block["idx"], "begin": begin}
        logging.info("vertical blocks: %s", vertical_blocks)

        #
        # We keep a set of all columns that are covered by a vertical block:
        # those columns will not be involved in the U[r,c] variables and in the
        # corresponding covering constraints
        logging.info("loop: covered by vertical blocks")
        covered_by_vertical_block = set()
        for begin, item in vertical_blocks.items():
            block = self.input_blocks[item['idx']]
            # logging.info("Vertical block: %s, block.str()")
            for col in range(block.start, block.end + 1):
                covered_by_vertical_block.add(col)
        logging.info(
            "Covered by vertical blocks: %s", covered_by_vertical_block)
        logging.info(
            "No. covered by vertical blocks: %s out of %s" % (len(covered_by_vertical_block), self.n_cols))

        # msa_positions is a list of all positions (r,c) that are required to be
        # covered. We exclude the positions covered by vertical blocks, since
        # they will be guaranteed to be covered, as an effect of the fact that
        # the corresponding C variables will be set to 1
        logging.debug("generating msa positions")
        msa_positions = [(r, c) for r in range(self.n_seqs)
                         for c in set(range(self.n_cols)) - covered_by_vertical_block]

        # We compute a dictionary called zones with key the column and value the
        # zone_id, that is a progressive id of the region, where each region is
        # a set of consecutive columns that are either disjoint or included in a vertical block.
        logging.info("zones")
        current_zone = 0
        zone = {0: 0}
        zone_boundaries = {0: {'start': 0, 'end': -1}}
        # We keep also a dictionary called zone_boundaries with key the zone_id
        # and values the first and last column of the zone
        for col in range(1, self.n_cols):
            if (col in covered_by_vertical_block) != (col - 1 in covered_by_vertical_block):
                # logging.debug("Boundary at %s: %s-%s", col, (col in covered_by_vertical_block),
                #               (col - 1 in covered_by_vertical_block))

                zone_boundaries[current_zone]['end'] = col - 1
                current_zone += 1
                zone_boundaries[current_zone] = {'start': col, 'end': -1}
            zone[col] = current_zone
        zone_boundaries[current_zone]['end'] = self.n_cols - 1

        # We compute the set unsplit_blocks, corresponding to the
        # blocks that are disjoint from or contained in vertical blocks.
        # The blocks that we will encode with a C variable in the
        # ILP correspond to the blocks that are disjoint from vertical blocks or
        # are vertical blocks themselves.
        logging.debug("zones\n%s", zone)
        logging.debug("boundaries:\n%s", zone_boundaries)
        # We start with all vertical blocks, since they will be encoded with a C variable
        private_blocks = {
            (self.input_blocks[vertical_blocks[k]['idx']].start,
             self.input_blocks[vertical_blocks[k]['idx']].end,
             frozenset(self.input_blocks[vertical_blocks[k]['idx']].K)): self.input_blocks[vertical_blocks[k]['idx']] for k in vertical_blocks}
        logging.debug("private vertical blocks: %s", private_blocks)

        for idx, block in enumerate(self.input_blocks):
            logging.debug("Analyzing %s out of %s. Size=%sx%s. Block=%s %s %s",
                         idx+1, n_blocks, len(block.K), block.end - block.start + 1, block.start, block.end, block.K)
            # Compute the zone of the boundaries of the block
            (zone_start, zone_end) = (zone[block.start], zone[block.end])
            logging.debug("Zones: %s-%s" % (zone_start, zone_end))
            # current_columns is the set of columns covered by the current block
            if (zone_start == zone_end):
                # Since the current block is contained in a single zone, it is
                # included in a vertical block or disjoint from vertical blocks.
                # In both cases, we encode it with a C variable
                logging.debug("unsplit block: %s %s %s" %
                             (block.start, block.end, block.K))
                # If the block is included in a vertical block, we can
                # discard it.
                # Notice that all maximal vertical blocks have been added to
                # private blocks before this for loop.
                if block.start not in covered_by_vertical_block:
                    private_blocks[(block.start, block.end,
                                    frozenset(block.K))] = block
                    logging.debug("Added %s -> %s", (block.start, block.end,
                                                    frozenset(block.K)), (block.start, block.end, block.K))
                # logging.debug("disjoint vertical block: %s", block.str())
            else:
                # The current block overlaps the vertical blocks.
                # We need to decompose it into parts that are disjoint from the
                # vertical blocks.
                logging.info("splitting block: %s %s %s" %
                             (block.start, block.end, block.K))

                # logging.debug("block %s is not disjoint", block.str())
                # logging.debug("zone_start: %s, zone_end: %s",
                #               zone_start, zone_end)
                # logging.debug("covered_by_vertical_block: %s",
                #               covered_by_vertical_block)

                # logging.debug(
                #     "block %s intersects with at least two vertical blocks", block.str())

                for zone_id in range(zone_start, zone_end + 1):
                    begin, end = max(zone_boundaries[zone_id]['start'],block.start), min(zone_boundaries[zone_id]['end'], block.end)
                    if begin not in covered_by_vertical_block:
                        # the current zone is not covered by a vertical block
                        if block.end < end:
                            # this is the last zone of the block, so the
                            # last position is the end of the block, not the
                            # end of the zone
                            end = block.end
                        label = str(self.msa[block.K[0], begin:end+1].seq)
                        new_block = Block(block.K, begin, end, label)
                        # check correct label of the block for all rows
                        for r in new_block.K:
                            for pos_block, c in enumerate(range(new_block.start, new_block.end + 1)):

                                # sanity check of labels
                                r = int(r)
                                char_msa = self.msa[r,c].upper()
                                char_block = new_block.label[pos_block].upper() 
                                if char_msa != char_block:
                                    logging.info("incorrect new_block covering position (%s,%s): %s" % (r,c,new_block.str()))
                                    logging.info("label of new_block should be %s", str(self.msa[r, new_block.start:new_block.end + 1].seq))


                        logging.debug("considering adding: %s %s %s" %
                                      (new_block.start, new_block.end, new_block.K))
                        private_blocks[(new_block.start, new_block.end,
                                        frozenset(new_block.K))] = new_block
                        logging.debug("Added %s -> %s", (new_block.start, new_block.end,
                                                        frozenset(new_block.K)),  (new_block.start, new_block.end, new_block.K))

                        # logging.debug("Adding private block: %s to %s" % (new_block.str(), private_blocks))
        logging.debug("Private blocks: %s", private_blocks)

        # Remove duplicates
        good_blocks = list(private_blocks.values())
        vertical_blocks = [idx for idx, block in enumerate(good_blocks) if (
            zone[block.start] == zone[block.end]) and block.start in covered_by_vertical_block]
        del private_blocks
        logging.info(
            f"Total number of distinct blocks: {len(good_blocks)}")
        logging.info("Vertical blocks: %s", len(vertical_blocks))

        logging.info("collecting c_variables")
        c_variables = list(range(len(good_blocks)))

        # covering_by_position is a dictionary with key (r,c) and value the list
        # of indices of the blocks that cover the position (r,c)
        logging.info("covering by position")
        n_cvars = len(c_variables)
        logging.info(f"MSA: {self.n_seqs} x {self.n_cols}")
        for idx in c_variables:
            block = good_blocks[idx]
            logging.debug("block: %s / %s" % (idx, n_cvars))
            logging.debug("Adding %s %s %s" %
                          (block.start, block.end, block.K))
            for r in block.K:
                for pos_block, c in enumerate(range(block.start, block.end + 1)):
                    covering_by_position[(r, c)].append(idx)

                    # sanity check of labels
                    r = int(r)
                    start = int(block.start)
                    end = int(block.end)
                    char_msa = self.msa[r,c].upper()
                    char_block = block.label[pos_block].upper() 
                    if char_msa != char_block:
                        logging.info("incorrect labeled block covering position (%s,%s): %s" % (r,c,block.str()))

        logging.info("Covering not vertical")
        for r in range(self.n_seqs):
            for c in set(range(self.n_cols)) - covered_by_vertical_block:
                logging.debug("Covering position: %s %s %s" %
                              (r, c, [(idx, good_blocks[idx].str()) for idx in covering_by_position[(r, c)]]))
                if len(covering_by_position[(r, c)]) == 0:
                    logging.debug("Uncovered position! %s %s" % (r, c))

        vertical_covered = set()
        for idx in vertical_blocks:
            logging.debug("Vertical block fixed: idx:%s %s" %
                          (idx, good_blocks[idx].str()))
            for col in range(good_blocks[idx].start, good_blocks[idx].end + 1):
                vertical_covered.add(col)
        logging.debug("Vertical covered: %s" %
                     covered_by_vertical_block.difference(vertical_covered))
        logging.debug("Vertical covered: %s" %
                     covered_by_vertical_block == vertical_covered)

        tf = time.time()
        times["init"] = round(tf - ti, 3)

        # Create the model
        ti = time.time()
        logging.info("initializing model pangeblocks")
        model = gp.Model("pangeblocks")

        # Threads
        model.setParam(GRB.Param.Threads, 16)

        # Time Limit
        model.setParam(GRB.Param.TimeLimit, float(60*self.time_limit))

        # define variables
        # C(b) = 1 if block b is selected
        # U(r,c) = 1 if position (r,c) is covered by at least one block
        # S(r,c) = 1 if position (r,c) is covered by a 1-cell block
        logging.info("adding C variables to the model")
        C = model.addVars(c_variables,
                          vtype=GRB.BINARY, name="C")
        for block in c_variables:
            logging.debug(
                "variable:C(%s) = %s" % (block, good_blocks[block].str()))

        logging.info("adding U variables to the model")
        U = model.addVars(msa_positions, vtype=GRB.BINARY, name="U")
        for pos in msa_positions:
            logging.debug("variable:U(%s)", pos)

        # Constraints
        # All C(b) variables corresponding to vertical blocks are set to 1
        logging.info("adding constraint: C=1 for vertical blocks ")
        for idx in vertical_blocks:
            model.addConstr(C[idx] == 1,
                            name=f"vertical_constraint_good_blocks({idx})")
            logging.debug("Vertical block fixed: %s %s" %
                         (idx, good_blocks[idx].str()))

        logging.info("adding constraints for each (r,c) position of the MSA")
        for r, c in msa_positions:
            logging.debug("MSA position (%s,%s)" % (r, c))
            blocks_rc = covering_by_position[(r, c)]
            if len(blocks_rc) > 0:
                # 1. U[r,c] = 1 implies that at least one block covers the position
                logging.debug("constraint 1")
                logging.debug("blocks_rc: %s" % blocks_rc)
                logging.debug("C variables: %s" % [C[i] for i in blocks_rc])
                model.addConstr(
                    U[r, c] <= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint1({r},{c})"
                )
                # 2. each position in the MSA is covered at most by one block
                logging.debug("constraint 2")
                model.addConstr(
                    1 >= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint2({r},{c})"
                )
                # logging.debug("constraint2(%s,%s) covered by %s" %
                #               (r, c, blocks_rc))
                # 3. each position of the MSA is covered AT LEAST by one block
                logging.debug("constraint 3")
                model.addConstr(U[r, c] >= 1, name=f"constraint3({r},{c})")
                # logging.debug("constraint3(%s,%s" % (r, c))
        tf = time.time()
        times["constraints1-2-3"] = round(tf - ti, 3)

        ti = time.time()
        # Objective function
        logging.info("setting objective function to %s", self.obj_function)
        if self.obj_function == "nodes":
            # # minimize the number of blocks (nodes)
            model = loss_nodes(model, vars=C)

        elif self.obj_function == "strings":
            # # minimize the total length of the graph (number of characters)
            model = loss_strings(model, vars=C, blocks=good_blocks, c_variables=c_variables)

        elif self.obj_function == "weighted":
            # # minimize the number of blocks penalizing shorter blocks
            model = loss_weighted(model, vars=C, blocks=good_blocks, c_variables=c_variables, 
                                  penalization=self.penalization, min_len=self.min_len)
            logging.info(f"penalization: {self.penalization}")
            logging.info(f"minimum length: {self.min_len}")

        elif self.obj_function == "depth":
            # # minimize the number blocks covering less than k sequences
            model = loss_depth(model, vars=C, blocks=good_blocks, c_variables=c_variables, 
                                  penalization=self.penalization, min_coverage=self.min_coverage, n_seqs=self.n_seqs)
            logging.info(f"penalization: {self.penalization}")
            logging.info(f"minimum coverage: {self.min_coverage}")
        
        tf = time.time()
        times["objective function"] = round(tf - ti, 3)

        ti = time.time()
        

        # ---  Optimization ---  
        model.optimize()
        logging.info("End ILP")
        tf = time.time()
        times["optimization"] = round(tf - ti, 3)

        if self.path_save_ilp:
            logging.info("saving ILP model")
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
                logging.debug(
                    f"Optimal Solution: {k}, {good_blocks[k].str()}")
                optimal_coverage.append(good_blocks[k])
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
