import gurobipy as gp
from gurobipy import GRB

def loss(model, vars, blocks, c_variables, penalization, min_coverage, n_seqs):

    MIN_COVERAGE= min_coverage # penalize blocks covering less than MIN_COVERAGE % of the sequences
    PENALIZATION = penalization # costly than others
    model.setObjective(
        gp.quicksum(
            (1 if len(blocks[idx].K)/n_seqs > MIN_COVERAGE  else PENALIZATION)*vars[idx] 
            for idx in c_variables
        )
    )
    return model