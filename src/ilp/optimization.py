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
        self.time_limit = kwargs.get("time_limit", 180)
        logging.getLogger().setLevel(log_level)

    def __call__(self, return_times: bool = False):
        "Solve ILP formulation"
        times = dict()
        ti = time.time()
        
        all_blocks = self.input_blocks
        n_blocks = len(all_blocks)
        logging.info("number of blocks %s", n_blocks)

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
        for idx, block in enumerate(all_blocks):
            # logging.debug("input block: %s", block.str())
            if len(block.K) == self.n_seqs:
                if not block.i in left_maximal_vertical_blocks or block.j > left_maximal_vertical_blocks[block.i]["end"]:
                    left_maximal_vertical_blocks[block.i] = {
                        "idx": idx, "end": block.j}
        
        logging.info("loop: vertical blocks")
        vertical_blocks = {}
        for begin, block in left_maximal_vertical_blocks.items():
            if not block["end"] in vertical_blocks or begin < vertical_blocks[block["end"]]["begin"]:
                vertical_blocks[block["end"]] = {
                    "idx": block["idx"], "begin": begin}
        #
        # We keep a set of all columns that are covered by a vertical block:
        # those columns will not be involved in the U[r,c] variables and in the
        # corresponding covering constraints
        logging.info("loop: covered by vertical blocks")
        covered_by_vertical_block = set()
        for begin, item in vertical_blocks.items():
            block = all_blocks[item['idx']]
            # logging.info("Vertical block: %s, block.str()")
            for col in range(block.i, block.j + 1):
                covered_by_vertical_block.add(col)
        logging.info(
            "Covered by vertical blocks: %s", covered_by_vertical_block)
        logging.info(
            "No. covered by vertical blocks: %s out of %s" % (len(covered_by_vertical_block), self.n_cols))

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

        # We compute the set disjoint_vertical, corresponding to the
        # blocks that are disjoint from vertical blocks.
        # The blocks that we will encode with a C variable in the
        # ILP correspond to the blocks that are disjoint from vertical blocks or
        # are vertical blocks themselves.
        # logging.debug("zones\n%s", zone)
        # logging.debug("boundaries:\n%s", zone_boundaries)
        intersecting_vertical_list = []
        disjoint_vertical = set()
        private_blocks = []
        for idx, block in enumerate(all_blocks):
            logging.info("analyzing %s out of %s", idx, n_blocks)
            # Compute the zone of the boundaries of the block
            (zone_start, zone_end) = (zone[block.i], zone[block.j])
            # current_columns is the set of columns covered by the current block
            if (zone_start == zone_end):
                # The current block is contained in a single zone or it is
                # included in a vertical block.
                # In the latter case, we do not need to encode it with a C
                # variable, so we discard the block
                if not zone_start in covered_by_vertical_block:
                    # The current block is disjoint from vertical blocks
                    logging.info("block disjoint from vertical blocks")
                    disjoint_vertical.add(idx)
                # logging.debug("disjoint vertical block: %s", block.str())
            else:
                # The current block overlaps the vertical blocks.
                # We need to check if the block intersects with at least two
                # vertical blocks: in that case we need to manually decompose
                # it.
                intersecting_vertical_list.append(idx)
                if zone_end - zone_start >= 3 or (zone_end - zone_start >= 2) and (zone_start in covered_by_vertical_block or zone_end in covered_by_vertical_block):
                    logging.info("block is not disjoint")
                    # logging.debug("block %s is not disjoint", block.str())
                    # logging.debug("zone_start: %s, zone_end: %s",
                    #               zone_start, zone_end)
                    # logging.debug("covered_by_vertical_block: %s",
                    #               covered_by_vertical_block)

                    # logging.debug(
                    #     "block %s intersects with at least two vertical blocks", block.str())
                    for zone_id in range(zone_start, zone_end + 1):
                        start, end = zone_boundaries[current_zone]['start'], zone_boundaries[current_zone]['end']
                        # LAST CHANGE
                        label = str(self.msa[block.K[0], begin:end+1].seq)
                        new_block = Block(block.K, begin, end, label)
                        private_blocks.append(new_block)
                        # logging.debug("Adding private block: %s to %s" % (new_block.str(), private_blocks))

        intersecting_vertical = set(intersecting_vertical_list)
        first_private_block = len(all_blocks)
        # all_blocks += enumerate(private_blocks, first_private_block)
        all_blocks += private_blocks

        logging.info("collecting c_variables")
        c_variables = list(disjoint_vertical) + [item["idx"]
                                                 for item in vertical_blocks.values()] + list(range(first_private_block, len(all_blocks)))
        # msa_positions is a list of all positions (r,c) that are required to be
        # covered. We exclude the positions covered by vertical blocks, since
        # they will be guaranteed to be covered, as an effect of the fact that
        # the corresponding C variables will be set to 1
        logging.info("generating msa positions")
        msa_positions = [(r, c) for r in range(self.n_seqs)
                         for c in set(range(self.n_cols)) - covered_by_vertical_block]

        # covering_by_position is a dictionary with key (r,c) and value the list
        # of indices of the blocks that cover the position (r,c)
        logging.info("covering by position")
        n_cvars=len(c_variables)
        for idx in c_variables:
            block = all_blocks[idx]
            logging.info("block: %s / %s" %  (idx,n_cvars))
            if idx not in intersecting_vertical and idx not in vertical_blocks:
                for r in block.K:
                    for c in range(block.i, block.j + 1):
                        covering_by_position[(r, c)].append(idx)

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
            logging.info(
                "variable:C(%s) = %s" % (block, all_blocks[block].str()))
        
        logging.info("adding U variables to the model")
        U = model.addVars(msa_positions, vtype=GRB.BINARY, name="U")
        for pos in msa_positions:
            logging.info("variable:U(%s)", pos)

        # Constraints
        # All C(b) variables corresponding to vertical blocks are set to 1
        logging.info("adding constraint: C=1 for vertical blocks ")
        for end, item in vertical_blocks.items():
            model.addConstr(C[item["idx"]] == 1,
                            name=f"vertical_constraint({item['idx']})")
        
        logging.info("adding constraints for each (r,c) position of the MSA")
        for r, c in msa_positions:
            blocks_rc = covering_by_position[(r, c)]
            if len(blocks_rc) > 0:
                # 1. U[r,c] = 1 implies that at least one block covers the position
                logging.info("constraint 1")
                model.addConstr(
                    U[r, c] <= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint1({r},{c})"
                )
                # 2. each position in the MSA is covered at most by one block
                logging.info("constraint 2")
                model.addConstr(
                    1 >= gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint2({r},{c})"
                )
                # logging.debug("constraint2(%s,%s) covered by %s" %
                #               (r, c, blocks_rc))
                # 3. each position of the MSA is covered AT LEAST by one block
                logging.info("constraint 3")
                model.addConstr(U[r, c] >= 1, name=f"constraint3({r},{c})")
                # logging.debug("constraint3(%s,%s" % (r, c))
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
        logging.info("setting objective function to %s", self.obj_function)
        if self.obj_function == "nodes":
            # minimize the number of blocks (nodes)
            model.setObjective(C.sum("*"), GRB.MINIMIZE)
            # vars=[var for var in model.getVars() if var.VarName.startswith("C")]
            # model.setObjective(LinExpr([1. for _ in vars], [model.getVarByName(name) for name in vars]), GRB.MINIMIZE)
        elif self.obj_function == "strings":
            # minimize the total length of the graph (number of characters)
            model.setObjective(
                gp.quicksum(
                    all_blocks[idx].len()*C[idx]
                    for idx in c_variables
                ),
                GRB.MINIMIZE
            )
        elif self.obj_function == "weighted":
            # minimize the number of blocks penalizing shorter blocks
            MIN_LEN = self.min_len  # penalize blocks with label less than MIN_LEN
            PENALIZATION = self.penalization  # costly than the other ones
            model.setObjective(
                gp.quicksum(
                    (PENALIZATION if all_blocks[idx].len(
                    ) < MIN_LEN else 1)*C[idx]
                    for idx in c_variables
                ),
                GRB.MINIMIZE
            )

        tf = time.time()
        times["objective function"] = round(tf - ti, 3)

        ti = time.time()
        if self.obj_function == "weighted":
            logging.info(f"penalization: {self.penalization}")
            logging.info(f"minimum length: {self.min_len}")

#        model.optimize()
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
