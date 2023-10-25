import gurobipy as gp
from gurobipy import GRB

def loss(model, vars, blocks, c_variables, msa, start_column):
    model.setObjective(
        gp.quicksum(
            len(msa[int(blocks[idx].K[0]), int(blocks[idx].start - start_column):int(blocks[idx].end + 1 - start_column)].seq.replace("-","")) * vars[idx]
            for idx in c_variables
        ),
        GRB.MINIMIZE
    )
    
    return model