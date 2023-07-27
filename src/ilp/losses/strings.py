import gurobipy as gp
from gurobipy import GRB
# FIXME: label from MSA
def loss(model, vars, blocks, c_variables):
    model.setObjective(
        gp.quicksum(
            blocks[idx].len() * vars[idx]
            for idx in c_variables
        ),
        GRB.MINIMIZE
    )
    
    return model