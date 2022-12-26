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

        # write idx for blocks 
        blocks = [(id_block, b.i, b.j) for id_block,b in enumerate(self.input_blocks)] # (K,i,j)
        block_by_id = {j: b for j,b in enumerate(blocks)}

        # blocks covering each position
        covering_by_position=defaultdict(list)
        for block in blocks:
            id_block,i,j = block
            K = [int(r) for r in id_block_to_K[id_block].split(",")]
            for r in K:
                
                for c in range(block[1],block[2]+1):
                    covering_by_position[(r,c)].append(id_block)

        # write idx for MSA positions (row, col)
        tf = time.time()
        times["init"] = round(tf - ti, 3)

        # Create the model
        ti = time.time()
        model = gp.Model("pangeblocks")

        # define variables
        for block in range(len(self.input_blocks)):
            logging.info(
                f"variable:C({block}) = {self.input_blocks[block].str()}")
        U = model.addVars(msa_positions, vtype=GRB.BINARY, name="U")
        for pos in msa_positions:
            logging.info(f"variable:U({pos})")

        # Constraints
        for r,c in tqdm(msa_positions):

            blocks_rc = covering_by_position[(r,c)] # TODO: save name of variables in covering by position 
            subset_C=[]
            for block_id in blocks_rc:
                b=block_by_id[block_id]
                subset_C.append(
                    C[b[0],b[1],b[2]]
                )
            
            if len(subset_C)>0:
                # import pdb; pdb.set_trace()
                ## 1. each position in the MSA is covered ONLY ONCE
                model.addConstr( U[r,c] <= gp.quicksum(subset_C), name=f"constraint1({r},{c})")
                
        tf = time.time()
        times["constraints1-2"] = round(tf - ti, 3)

        tf=time.time()
        times["constraints1-2"] = round(tf-ti,3)

        ## 3. overlapping blocks cannot be chosen
        # sort all blocks,
        ti = time.time()

        intersection = self._list_inter_blocks(self.input_blocks)

        for pos1,pos2 in tqdm(intersection):
            # if the blocks intersect, then create the restriction 
            block1=blocks[pos1]
            block2=blocks[pos2]
            K1,i1,j1=block1
            K2,i2,j2=block2
            name_constraint=f"constraint3({K1},{i1},{j1})-({K2},{i2},{j2})"
            model.addConstr(C[block1] + C[block2] <= 1 , name=name_constraint)

        tf = time.time()
        times["constraint3"] = round(tf - ti, 3)

        ti = time.time()
        # Objective function
        model.setObjective(C.sum('*','*','*'), GRB.MINIMIZE)

        model.optimize()
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
        for k,v in solution_C.items(): 
            id_block,i,j=k 
            if v > 0:
                K = (int(seq) for seq in id_block_to_K[id_block].split(","))
                label = id_block_to_labels[id_block]
                optimal_coverage.append(
                    Block(K,i,j,label)
                )
        tf=time.time()
        times["solution as blocks"] = round(tf-ti,3)

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
        blocks = [(idx,block) for idx,block in sorted(enumerate(list_blocks), key=lambda block: block[1].i)]
        
        # save pairs of indexes for the sorted blocks that intersect
        # (pos1,pos2): positions in the sorted list
        # (orig_pos1, orig_pos2): positions in the order of the input list
        intersections = [] 
        for pos1, pos_block1 in enumerate(blocks[:-1]):
            idx1, block1 = pos_block1
            orig_pos1 = idx1 
                        
            # compare against the next blocks in the sorted list 
            for rel_pos, pos_block2 in enumerate(blocks[pos1+1:]):
                idx2, block2 = pos_block2
                orig_pos2 = idx2
                pos2 = rel_pos + pos1 + 1
                
                if block1.j < block2.i: break # no intersection is possible
                
                # check for not empty intersection
                common_rows = list(set(block1.K).intersection(set(block2.K))) # intersection set K
                common_cols = self.get_common_cols(block1,block2) # intersection columns [i,j]

                if (common_rows and common_cols):
                    intersections.append((orig_pos1,orig_pos2))

        return intersections
