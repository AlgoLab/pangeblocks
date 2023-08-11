from .. import Block


import logging
logging.basicConfig(level=logging.INFO,
                    format='[Block Decomposition Row Maximal] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

# logging.Formatter(
#     fmt='%(asctime)s.%(msecs)03d',
#     datefmt='%Y-%m-%d,%H:%M:%S'
# )
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
    b1,b2=sorted([block1,block2], key=lambda b: (b.start,b.end))

    nb = [] # new blocks
    # not empty intersection
    common_rows = list(set(b1.K).intersection(set(b2.K)))
    common_cols = list(set(range(b1.start,b1.end+1)).intersection(set(range(b2.start,b2.end+1))))
    
    if common_rows and common_cols: 
        K = list(common_rows)
        K.sort()

        # Condition 1
        if b1.start == b2.start and b1.end < b2.end:
            logging.info("decomposition 1")
            # new blocks
            nb1 = b1#Block(b1.K, b1.start, b1.end, b1.label)
            nb2 = Block(b2.K, b1.end+1, b2.end)
            nb.extend([nb1, nb2])

        # Condition 2
        elif b1.start < b2.start and b2.end < b1.end:
            logging.info("decomposition 2")
            nb1 = Block(b1.K, b1.start, b2.start-1)
            nb2 = b2
            nb3 = Block(b1.K, b2.end+1, b1.end)
            nb.extend([nb1, nb2, nb3])

        # Condition 3
        elif b1.start < b2.start and b1.end == b2.end:
            logging.info("decomposition 3")
            nb1 = Block(b1.K, b1.start, b2.start-1)
            nb2 = b2
            nb.extend([nb1, nb2])

        # Condition4
        elif b1.start < b2.start and b2.start < b1.end and b1.end < b2.end:
            logging.info("decomposition 4")
            # option 1 
            nb1 = Block(b1.K, b1.start, b2.start-1)
            nb2 = b2
            nb.extend([nb1, nb2])

            # option 2
            nb1 = b1
            nb2 = Block(b2.K, b1.end+1, b2.end)
            nb.extend([nb1, nb2])
            
    return nb