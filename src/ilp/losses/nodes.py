import gurobipy as gp
from gurobipy import GRB

def loss(model, vars,):
    model.setObjective(vars.sum("*"), GRB.MINIMIZE)
    return model