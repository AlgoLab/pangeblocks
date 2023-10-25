import gurobipy as gp
from gurobipy import GRB

def loss(model, vars, blocks, c_variables, start_column, msa):
    # Given a block l=(K,b,e), the cost is w(l) = f(l)/|K|, where f(l)= (b-e+1) - #indels
    model.setObjective(
        gp.quicksum(
            len(msa[int(blocks[idx].K[0]), int(blocks[idx].start - start_column):int(blocks[idx].end + 1 - start_column)].seq.replace("-","")) / len(blocks[idx].K) *vars[idx] \
            for idx in c_variables
        )
    )
    return model