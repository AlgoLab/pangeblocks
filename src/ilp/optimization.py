"""Solve the ILP formulation:
Given an MSA, and a set of blocks such that the union of them covers all the MSA,
Find a non-overlapping set the blocks covering the entire MSA
# """
# from logging import Logger
# log = Logger(name="opt", level="DEBUG")

import gurobipy as gp
from gurobipy import GRB
from Bio import AlignIO
from pathlib import Path
from ..blocks import Block
from tqdm import tqdm 

class Optimization:
    
    def __init__(self, blocks, path_msa, path_save_ilp=None):

        self.blocks = blocks
        msa, n_seqs, n_cols = self.load_msa(path_msa)
        self.msa = msa
        self.n_seqs = n_seqs
        self.n_cols = n_cols
        self.path_save_ilp = path_save_ilp

    def __call__(self,):
        """Solve ILP formulation"""
        id_block_to_K = {id_block: ",".join([str(r) for r in block.K]) for id_block, block in enumerate(self.blocks) }
        id_block_to_labels = {id_block: b.label for id_block, b in enumerate(self.blocks)}

        # write idx for blocks 
        blocks = [(id_block, b.i, b.j) for id_block,b in enumerate(self.blocks)] # (K,i,j)

        # write idx for MSA positions (row, col)
        msa_positions = [(r,c) for r in range(self.n_seqs) for c in range(self.n_cols)] 

        # Create the model
        model = gp.Model("pangeblocks")

        # define variables
        C = model.addVars(blocks, vtype=GRB.BINARY, name="C")
        U = model.addVars(msa_positions, vtype=GRB.BINARY, name="U")

        # Constraints:
        for r,c in tqdm(msa_positions):

            # subset of blocks that covers the position [r,c]
            subset_C = [ C[id_block,i,j] for id_block,i,j in blocks if str(r) in id_block_to_K[id_block].split(",") and i<=c<=j ]

            if len(subset_C)>0:                
                ## 1. each position in the MSA is covered ONLY ONCE
                model.addConstr( U[r,c] <= sum(subset_C), name=f"constraint1({r},{c})")
                
                ## 2. each position of the MSA is covered AT LEAST by one block
                model.addConstr( U[r,c] >= 1, name=f"constraint2({r},{c})")

        ## 3. overlapping blocks cannot be chosen
        # sort all blocks, 
        blocks = sorted(blocks, key=lambda b: b[1]) # sort blocks by the starting position (K,start,end)

        # and analyze the intersections while update the constraints
        names_constraint3=[]
        for pos1,block1 in enumerate(blocks[:-1]):
            # compare against the next blocks in the sorted list
            for rel_pos, block2 in enumerate(blocks[pos1+1:]):
                pos2 = rel_pos + pos1 + 1
                if block1[2] < block2[1]: 
                    break # no intersection is possible 
                block2 = blocks[pos2]
                
                # check for not empty intersection, otherwise, skip to the next block  
                # note: set K is a string with the rows concatenated by a "," (due to Gurobi requirements to index the variables)
                id_block1 = block1[0]
                id_block2 = block2[0]
                block1_K = id_block_to_K[id_block1].split(",")
                block2_K = id_block_to_K[id_block2].split(",")

                # check for not empty intersection, otherwise skip to the next block1 in the list
                common_rows = list(set(block1_K).intersection(set(block2_K))) # intersection set K
                # common_cols = list(set(range(block1[1],block1[2]+1)).intersection(set(range(block2[1],block2[2]+1)))) # intersection columns [i,j]
                common_cols = self.get_common_cols(block1,block2)

                if (common_rows and common_cols):
                    
                    # if the blocks intersect, then create the restriction 
                    K1,i1,j1=block1
                    K2,i2,j2=block2
                    name_constraint=f"constraint3({K1},{i1},{j1})-({K2},{i2},{j2})"
                    model.addConstr(C[block1] + C[block2] <= 1 , name=name_constraint)
                    names_constraint3.append(name_constraint)

        # Objective function
        model.setObjective(C.sum('*','*','*'), GRB.MINIMIZE)

        model.optimize()

        if self.path_save_ilp: 
            Path(self.path_save_ilp).parent.mkdir(exist_ok=True, parents=True)
            model.write(self.path_save_ilp)

        try:                
            solution_C = model.getAttr("X", C)
        except:
            raise("No solution")
        
        # filter optimal coverage of blocks for the MSA
        optimal_coverage = []
        for k,v in solution_C.items(): 
            id_block,i,j=k 
            if v > 0:
                K = (int(seq) for seq in id_block_to_K[id_block].split(","))
                label = id_block_to_labels[id_block]
                optimal_coverage.append(
                    Block(K,i,j,label)
                )

        return optimal_coverage


    def load_msa(self, path_msa):
        "return alignment, number of sequences and columns"
        # load MSA
        align=AlignIO.read(path_msa, "fasta")
        n_cols = align.get_alignment_length()
        n_seqs = len(align)

        return align, n_seqs, n_cols

    def get_common_cols(self, block1,block2):
        if block1[2] < block2[1]:
            return []
        intervals = [(block1[1],block1[2]),(block2[1],block2[2])]
        start, end = intervals.pop()
        while intervals:
            start_temp, end_temp = intervals.pop()
            start = max(start, start_temp)
            end = min(end, end_temp)
        return [start, end]