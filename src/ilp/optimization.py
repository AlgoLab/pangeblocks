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

        # TODO: input_blocks must be the complete list of blocks that the ILP will use
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

        # TODO: Modify until here
        
        logging.info("collecting c_variables")
        c_variables = list(range(len(self.input_blocks))) # TODO: self.input_blocks must be the input to the ILP

        # covering_by_position is a dictionary with key (r,c) and value the list
        # of indices of the blocks that cover the position (r,c)
        logging.info("covering by position")
        n_cvars = len(c_variables)
        logging.info(f"MSA: {self.n_seqs} x {self.n_cols}")
        for idx in c_variables:
            block = self.input_blocks[idx]
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

        # logging.info("Covering not vertical")
        # for r in range(self.n_seqs):
        #     for c in set(range(self.n_cols)) - covered_by_vertical_block:
        #         logging.debug("Covering position: %s %s %s" %
        #                       (r, c, [(idx, self.input_blocks[idx].str()) for idx in covering_by_position[(r, c)]]))
        #         if len(covering_by_position[(r, c)]) == 0:
        #             logging.debug("Uncovered position! %s %s" % (r, c))

        # vertical_covered = set()
        # for idx in vertical_blocks:
        #     logging.debug("Vertical block fixed: idx:%s %s" %
        #                   (idx, self.input_blocks[idx].str()))
        #     for col in range(self.input_blocks[idx].start, self.input_blocks[idx].end + 1):
        #         vertical_covered.add(col)
        # logging.debug("Vertical covered: %s" %
        #              covered_by_vertical_block.difference(vertical_covered))
        # logging.debug("Vertical covered: %s" %
        #              covered_by_vertical_block == vertical_covered)

        tf = time.time()
        times["init"] = round(tf - ti, 3)



        #  ------------------------
        #  --- Create the model ---
        #  ------------------------
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
                "variable:C(%s) = %s" % (block, self.input_blocks[block].str()))

        logging.info("adding U variables to the model")
        U = model.addVars(msa_positions, vtype=GRB.BINARY, name="U")
        for pos in msa_positions:
            logging.debug("variable:U(%s)", pos)

        #  ------------------------
        #  Constraints
        #  ------------------------
        # All C(b) variables corresponding to vertical blocks are set to 1
        logging.info("adding constraint: C=1 for vertical blocks ")
        for idx in vertical_blocks:
            model.addConstr(C[idx] == 1,
                            name=f"vertical_constraint_good_blocks({idx})")
            logging.debug("Vertical block fixed: %s %s" %
                         (idx, self.input_blocks[idx].str()))

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

        #  ------------------------
        #  Objective function
        #  ------------------------
        logging.info("setting objective function to %s", self.obj_function)
        if self.obj_function == "nodes":
            # # minimize the number of blocks (nodes)
            model = loss_nodes(model, vars=C)

        elif self.obj_function == "strings":
            # # minimize the total length of the graph (number of characters)
            model = loss_strings(model, vars=C, blocks=self.input_blocks, c_variables=c_variables)

        elif self.obj_function == "weighted":
            # # minimize the number of blocks penalizing shorter blocks
            model = loss_weighted(model, vars=C, blocks=self.input_blocks, c_variables=c_variables, 
                                  penalization=self.penalization, min_len=self.min_len)
            logging.info(f"penalization: {self.penalization}")
            logging.info(f"minimum length: {self.min_len}")

        elif self.obj_function == "depth":
            # # minimize the number blocks covering less than k sequences
            model = loss_depth(model, vars=C, blocks=self.input_blocks, c_variables=c_variables, 
                                  penalization=self.penalization, min_coverage=self.min_coverage, n_seqs=self.n_seqs)
            logging.info(f"penalization: {self.penalization}")
            logging.info(f"minimum coverage: {self.min_coverage}")
        
        tf = time.time()
        times["objective function"] = round(tf - ti, 3)

        
        #  ------------------------
        #  Solve or save the ILP model
        #  ------------------------
        
        # TODO: use flag to store the ILP instead of solving it
        ti = time.time()
        if self.path_save_ilp:
            logging.info("saving ILP model")
            Path(self.path_save_ilp).parent.mkdir(exist_ok=True, parents=True)
            model.write(self.path_save_ilp)
            tf = time.time()
            times["save-ilp"] = round(tf - ti, 3)
        else:
            logging.info("solving ILP model")
            model.optimize()
            logging.info("End ILP")
            tf = time.time()
            times["optimization"] = round(tf - ti, 3)

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
