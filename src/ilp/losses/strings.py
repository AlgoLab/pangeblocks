import gurobipy as gp
from gurobipy import GRB

def loss(model, vars, blocks, c_variables):
    model.setObjective(
        gp.quicksum(
            len(blocks[idx].label.replace("-","")) * vars[idx]
            for idx in c_variables
        ),
        GRB.MINIMIZE
    )
    
    return model