"""Solve the ILP formulation:
Given an MSA, and a set of blocks such that the union of them covers all the MSA,
Find a non-overlapping set the blocks covering the entire MSA
# """
import os
import sys
from sys import getsizeof, stderr
from collections import defaultdict
from blocks import Block
from pathlib import Path
from Bio import AlignIO
from gurobipy import GRB, LinExpr
import gurobipy as gp
from .losses import (
    loss_nodes,
    loss_strings,
    loss_weighted,
    loss_depth
)

import logging


def load_submsa(filename, start_column=0, end_column=-1):
    "Return MSA from start_column to end_column (both included)"
    # load MSA
    msa = AlignIO.read(filename, "fasta")
    
    # filter sub-MSA if start/end columns are given
    if start_column>0 or end_column!=-1:
        # get last column
        # n_seqs = len(align)
        n_cols = msa.get_alignment_length()
        assert start_column < n_cols and end_column < n_cols, f"start_column={start_column}, end_column={end_column}. Must be < {n_cols} (number of columns in the MSA)"
        msa = msa[:, start_column:end_column+1] # end_column included
    # logging.info(f"msa ({start_column}, {end_column}): {msa}")
    return msa

class Optimization:
    "Generates/solves ILP model for a subMSA, from column start_column to end_column"
    def __init__(self, blocks: list, path_msa: str, 
                 start_column: int, end_column: int,
                 path_save_ilp: str, log_level=logging.INFO, **kwargs):
        
        self.input_blocks=blocks
        self.path_msa=path_msa
        self.start_column=start_column
        self.end_column=end_column
        self.path_save_ilp = path_save_ilp
        logging.basicConfig(level=log_level,
                    format='[Optimization] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')
        # logging.getLogger().setLevel(log_level)

        # Each block is a tuple (K, i, j, label) where
        # K is a tuple of rows,
        # i and j are the first and last column
        # label is the string of the block
        msa = load_submsa(path_msa, start_column=start_column, end_column=end_column)
        self.msa = msa
        self.n_seqs = len(msa)
        self.n_cols = msa.get_alignment_length()
        
        # ILP params       
        self.obj_function = kwargs.get("obj_function", "nodes")
        self.penalization = kwargs.get("penalization", 1)
        self.min_len = kwargs.get("min_len", 1)
        self.min_coverage = kwargs.get("min_coverage", 1)
        self.time_limit = kwargs.get("time_limit", 180)
        self.threads_ilp = kwargs.get("threads_ilp", 4)
        
    def __call__(self, solve_ilp: bool = False):
        """Generates the ILP model, and solves it if return_model is False

        Args:
            return_times (bool, optional): return detailed time of the execution. Defaults to False.
            solve_ilp (bool, optional): try to find ptimal solution if True, otherwise just save the ILP formulation. Defaults to False.

        Returns: 
            list with blocks in the optimal solution
        """        
        n_blocks = len(self.input_blocks)
        logging.info(f"Number of blocks ilp {n_blocks} ({self.start_column},{self.end_column})")
        for idx, block in enumerate(self.input_blocks):
            logging.debug("input block %s: %s" % (idx, block.str()))

        # covering_by_position is a dictionary with key (r,c) and value the list
        # of indices of the blocks that include the position (r,c)
        # To speed up the process, we iterate over the blocks and append the
        # block to all the positions it covers
        covering_by_position = defaultdict(list)
        
        logging.info(f"collecting c_variables ({self.start_column},{self.end_column})")
        c_variables = list(range(len(self.input_blocks)))
        logging.info(f"Size bytes of c_variables list {sys.getsizeof(c_variables)} ({self.start_column},{self.end_column})")

        # covering_by_position is a dictionary with key (r,c) and value the list
        # of indices of the blocks that cover the position (r,c)
        logging.info("covering by position")
        n_cvars = len(c_variables)
        logging.info(f"Number of C variables {n_cvars} ({self.start_column},{self.end_column})")
        logging.info(f"MSA: {self.n_seqs} x {self.n_cols} ({self.start_column},{self.end_column})")
        for idx in c_variables:
            block = self.input_blocks[idx]
            logging.debug("block: %s / %s" % (idx, n_cvars))
            logging.debug("Adding %s %s %s" %
                          (block.start, block.end, block.K))
            for r in block.K:
                for pos_block, c in enumerate(range(block.start, block.end + 1)):
                    covering_by_position[(r, c)].append(idx)

                    # # sanity check of labels
                    # r = int(r)
                    # start = int(block.start)
                    # end = int(block.end)
                    # char_msa = self.msa[r,c].upper()
                    # char_block = block.label[pos_block].upper() 
                    # if char_msa != char_block:
                    #     logging.info("incorrect labeled block covering position (%s,%s): %s" % (r,c))
        logging.info(f"Size bytes of covering_by_position dict {sys.getsizeof(covering_by_position)} ({self.start_column},{self.end_column})")

        #  ------------------------
        #  --- Create the model ---
        #  ------------------------
        logging.info(f"initializing model pangeblocks ({self.start_column},{self.end_column})")
        model = gp.Model("pangeblocks")
        logging.info(f"Size bytes of model initialized {sys.getsizeof(model)} ({self.start_column},{self.end_column})")
        # Threads
        model.setParam(GRB.Param.Threads, self.threads_ilp)
        
        # Time Limit
        model.setParam(GRB.Param.TimeLimit, float(60*self.time_limit))

        # define variables
        # C(b) = 1 if block b is selected
        logging.info(f"adding C variables to the model ({self.start_column},{self.end_column})")
        C = model.addVars(c_variables,
                          vtype=GRB.BINARY, name="C")
        logging.info(f"added C variables to the model ({self.start_column},{self.end_column})")
        logging.info(f"Size bytes of C variables {sys.getsizeof(C)} ({self.start_column},{self.end_column})")        

        for block in c_variables:
            logging.debug(
                "variable:C(%s) = %s" % (block, self.input_blocks[block]))

        # logging.info(f"adding U variables to the model ({self.start_column},{self.end_column})")
        # FIXME: exclude positions covered by missing blocks
        msa_positions = [(r,c + self.start_column) for r in range(self.n_seqs) for c in range(self.n_cols)]
    
        #  ------------------------
        #  Constraints
        #  ------------------------
        logging.info(f"adding constraints for each (r,c) position of the MSA ({self.start_column},{self.end_column})")
        for r, c in msa_positions:
            logging.debug("MSA position (%s,%s)" % (r, c))
            blocks_rc = covering_by_position[(r, c)]
            if len(blocks_rc) > 0:
                logging.debug(f"position [{r},{c}] is covered by {len(blocks_rc)} blocks ({self.start_column},{self.end_column})")
        
                # 1. each position in the MSA is covered by one block
                logging.debug(f"start constraint position [{r},{c}] ({self.start_column},{self.end_column})")
                model.addConstr(
                    1 == gp.quicksum([C[i] for i in blocks_rc]), name=f"constraint({r},{c})"
                )
                logging.debug(f"end constraint position [{r},{c}] ({self.start_column},{self.end_column})")
        logging.info(f"added constraints for each (r,c) position of the MSA ({self.start_column},{self.end_column})")
        
        model.update()  # to make the constraint take effect https://www.gurobi.com/documentation/current/refman/py_constr.html

        # Get a list of all constraints in the model
        all_constraints = model.getConstrs()

        # Calculate the size of constraints in bytes
        total_size_bytes = 0
        for constr in all_constraints:
            # Use the __sizeof__() method to get the size of the constraint in bytes
            size_bytes = getsizeof(constr) #.__sizeof__()
            total_size_bytes += size_bytes

        # size_constraints=sum([constr.__sizeof__() for constr in model.getConstrs()])
        logging.info(f"Size bytes of constraints {total_size_bytes} ({self.start_column},{self.end_column})")
        
        #  ------------------------
        #  Objective function
        #  ------------------------
        logging.info("##############################################################")
        logging.info(f"set objective function ({self.start_column},{self.end_column})")
        logging.info("setting objective function to %s", self.obj_function)
        if self.obj_function == "nodes":
            # # minimize the number of blocks (nodes)
            model = loss_nodes(model, vars=C)

        elif self.obj_function == "strings":
            # # minimize the total length of the graph (number of characters)
            for idx in c_variables:
                logging.debug(f"{self.msa}")
                logging.debug(f"block: {self.input_blocks[idx]}")
                k0=int(self.input_blocks[idx].K[0])
                logging.debug(f"K[0]: {k0}")
                sc = int(self.input_blocks[idx].start - self.start_column)
                logging.debug(f"start column: {sc}")
                ec = int(self.input_blocks[idx].end+1-self.start_column)
                logging.debug(f"end column: {ec}")
                logging.debug(f"msa[{k0},{sc}:{ec}]")
                print(self.msa[k0])
                logging.debug(f"type {type(self.msa[int(k0)])}")
                logging.debug(f"type {type(self.msa[int(k0),sc:ec])}")

                logging.debug(f"msa[{k0},{sc}:{ec}]: {self.msa[k0,sc:ec]}")
                # logging.info(f"{self.msa[self.input_blocks[idx].K[0], self.input_blocks[idx].start - self.start_column:self.input_blocks[idx].end+1-self.start_column]}")
                # logging.info(f"{self.msa[self.input_blocks[idx].K[0], self.input_blocks[idx].start - self.start_column:self.input_blocks[idx].end+1-self.start_column]}")
            
            model = loss_strings(model, vars=C, blocks=self.input_blocks, c_variables=c_variables, msa=self.msa, start_column=self.start_column)

        elif self.obj_function == "weighted":
            # # minimize the number of blocks penalizing shorter blocks
            model = loss_weighted(model, vars=C, blocks=self.input_blocks, c_variables=c_variables, 
                                  penalization=self.penalization, min_len=self.min_len, msa=self.msa, start_column=self.start_column)
            logging.info(f"penalization: {self.penalization}")
            logging.info(f"minimum length: {self.min_len}")

        elif self.obj_function == "depth":
            # # minimize the number blocks covering less than k sequences
            model = loss_depth(model, vars=C, blocks=self.input_blocks, c_variables=c_variables, 
                                  penalization=self.penalization, min_coverage=self.min_coverage, n_seqs=self.n_seqs)
            logging.info(f"penalization: {self.penalization}")
            logging.info(f"minimum coverage: {self.min_coverage}")
        logging.info(f"setted objective function ({self.start_column},{self.end_column})")
        #  ------------------------
        #  Solve or save the ILP model
        #  ------------------------
        if solve_ilp:
            logging.info(f"Start ILP ({self.start_column},{self.end_column})")
            model.optimize()
            logging.info(f"End ILP ({self.start_column},{self.end_column})")
            
            try:
                solution_C = model.getAttr("X", C)
            except:
                raise ("No solution")

            # filter optimal coverage of blocks for the MSA
            optimal_coverage = []
            for k, v in solution_C.items():

                if v > 0:
                    logging.debug(
                        f"Optimal Solution: {k}, {self.input_blocks[k]}")
                    optimal_coverage.append(self.input_blocks[k])
            
        if self.path_save_ilp:
            logging.info("saving ILP model")
            Path(self.path_save_ilp).parent.mkdir(exist_ok=True, parents=True)
            model.write(self.path_save_ilp)

        return optimal_coverage
        