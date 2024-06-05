from .. import Block

import logging
logging.basicConfig(level=logging.INFO,
                    format='[Complete Decomposition] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

def block_decomposition(block1: Block, block2: Block):
    """Decompose 2 blocks based on their intersection

    Args:
        block1 (Block): a block
        block2 (Block): another block

    Returns:
        list: blocks decomposed from the intersection. If input blocks
                does not intersect, the output is an empty list.
    """
    # sort blocks (left most first)
    b1,b2=sorted([block1,block2], key=lambda b: b.start)

    blocks=[]
    # not empty intersection
    common_rows = list(set(b1.K).intersection(set(b2.K)))
    common_cols = list(set(range(b1.start,b1.end+1)).intersection(set(range(b2.start,b2.end+1))))

    if common_rows and common_cols:
        K = list(common_rows)
        K.sort()

        if b1.start < b2.start:                                   # manuscript eq (1)
            blocks.append( Block(b1.K, b1.start, b2.start - 1) )     

        if b1.end < b2.end:
            blocks.append( Block(b2.K, b1.end + 1 , b2.end) )     # manuscript eq (2)
        
        if b1.end > b2.end:
            blocks.append( Block(b1.K, b2.end + 1 , b1.end) )     # manuscript eq (3)

        # extra blocks, part of the complete decomposition 
        K_inter = set(b1.K).intersection(b2.K) 
        K1_minus_K2 = set(b1.K) - set(b2.K)
        K2_minus_K1 = set(b2.K) - set(b1.K)

        if K_inter:                                              # manuscript eq (4)
            blocks.append( Block(K_inter, b1.start, b1.end) )  

        if K1_minus_K2:                                          # manuscript eq (5)
            blocks.extend([
                Block(K1_minus_K2, b1.start, b1.end),
                Block(K1_minus_K2, b2.start, b2.end)
            ])
        
        if K2_minus_K1:                                          # manuscript eq (6)
            blocks.extend([
                Block(K2_minus_K1, b1.start, b1.end),
                Block(K2_minus_K1, b2.start, b2.end)
            ])
        
        if K_inter and b2.start < b1.end:                       # manuscript eq (7)
            blocks.append(
                Block(K_inter, b2.start, b1.end)
            )

        if K1_minus_K2 and b2.start < b1.end:                   # manuscript eq (8) 
            blocks.append(
                Block(K1_minus_K2, b2.start, b1.end)
            )
            
        if K2_minus_K1 and b2.start < b1.end:                  # manuscript eq (9) 
            blocks.append(
                Block(K2_minus_K1, b2.start, b1.end)
            )
 
    return blocks