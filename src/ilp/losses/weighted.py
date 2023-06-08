import gurobipy as gp
from gurobipy import GRB

def loss(model, vars, blocks, c_variables, penalization, min_len):
    PENALIZATION = penalization
    MIN_LEN = min_len
    
    model.setObjective(
                gp.quicksum(
                    (PENALIZATION if len(blocks[idx].label.replace("-","")) < MIN_LEN else 1)*vars[idx]
                    for idx in c_variables
                ),
                GRB.MINIMIZE
    )
    return model