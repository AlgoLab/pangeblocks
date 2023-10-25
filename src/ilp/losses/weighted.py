import gurobipy as gp
from gurobipy import GRB

def loss(model, vars, blocks, c_variables, penalization, min_len, msa, start_column):
    PENALIZATION = penalization
    MIN_LEN = min_len
    
    model.setObjective(
                gp.quicksum(
                    (PENALIZATION if len(msa[int(blocks[idx].K[0]), int(blocks[idx].start - start_column):int(blocks[idx].end + 1 - start_column)].seq.replace("-","")) <= MIN_LEN else 1)*vars[idx]
                    for idx in c_variables
                ),
                GRB.MINIMIZE
    )
    return model