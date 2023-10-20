import gurobipy as gp
from gurobipy import GRB

def loss(model, vars, blocks, c_variables, start_column, msa, n_seqs):

    model.setObjective(
        gp.quicksum(
            len(msa[int(blocks[idx].K[0]), int(blocks[idx].start - start_column):int(blocks[idx].end + 1 - start_column)].seq.replace("-","")) * \
                len(blocks[idx].K)/n_seqs *vars[idx] \
            for idx in c_variables
        )
    )
    return model