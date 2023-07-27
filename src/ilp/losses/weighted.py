import gurobipy as gp
from gurobipy import GRB
# FIXME: label from MSA
def loss(model, vars, blocks, c_variables, penalization, min_len):
    PENALIZATION = penalization
    MIN_LEN = min_len
    
    model.setObjective(
                gp.quicksum(
                    (PENALIZATION if blocks[idx].len() < MIN_LEN else 1)*vars[idx]
                    for idx in c_variables
                ),
                GRB.MINIMIZE
    )
    return model