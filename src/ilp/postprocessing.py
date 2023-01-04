"""Postprocessing of optimal solution
Glue blocks that are consecutive and shares the same rows

- In the optimal solution, two blocks sharing the same set of sequences
are either consecutive or disjoint
"""
from collections import defaultdict
from ..blocks import Block

def postprocessing(list_blocks: list[Block]) -> list[Block]: 

    new_blocks = []
    blocks_by_setK = defaultdict(list)
    
    for idx, block in enumerate(list_blocks):
        blocks_by_setK[block.K].append(idx)

    # glue consecutive blocks
    for setK, idx_blocks in blocks_by_setK.items():
    
        blocks = [(idx,list_blocks[idx]) for idx in idx_blocks]
        blocks = sorted(blocks, key=lambda b: b[1].i)

        consecutive_blocks = [] # temporary list to collect consecutive blocks
        for _, block in blocks: 
            print(block)
            if consecutive_blocks: # not empty
                last_block = consecutive_blocks[-1]

                # check if new block is consecutive to the last block in the list
                # blocks are sorted by starting position
                if last_block.j == block.i - 1:
                    consecutive_blocks.append(block)
                else: 
                    start = consecutive_blocks[0].i
                    end   = consecutive_blocks[-1].j
                    label = "".join([b.label for b in consecutive_blocks])
                    new_blocks.append(
                        Block(setK, start, end, label)
                    )
            else:
                consecutive_blocks.append(block)

        # after last block
        start = consecutive_blocks[0].i
        end   = consecutive_blocks[-1].j
        label = "".join([b.label for b in consecutive_blocks])
        new_blocks.append(
            Block(setK, start, end, label)
        )

    return new_blocks