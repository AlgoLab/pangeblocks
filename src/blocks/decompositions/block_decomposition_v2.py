from .block import Block

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
    b1,b2=sorted([block1,block2], key=lambda b: (b.i,b.j))

    nb = [] # new blocks
    # not empty intersection
    common_rows = list(set(b1.K).intersection(set(b2.K)))
    common_cols = list(set(range(b1.i,b1.j+1)).intersection(set(range(b2.i,b2.j+1))))

    if common_rows and common_cols:
        K = list(common_rows)
        K.sort()

        # Condition 1 
        if b1.i == b2.i and b1.j < b2.j:
            # new blocks
            nb1 = Block(b2.K, b1.i, b1.j, b1.label)
            nb2 = Block(b2.K, b1.j+1, b2.j, b2.label[b1.j-b1.i+1:])
            nb3 = Block(set(b1.K)-set(b2.K), b1.i, b1.j, b1.label) 
            nb.extend([nb1, nb2, nb3])

        # Condition 2
        elif b1.i < b2.i and b2.j < b1.j:
            nb1 = Block(b1.K, b1.i, b2.i-1, b1.label[:b2.i-1-b1.i+1])
            nb2 = Block(b1.K, b2.i, b2.j, b2.label)
            nb3 = Block(b1.K, b2.j+1, b1.j, b1.label[b2.j-b1.i+1:])
            nb4 = Block(set(b2.K)-set(b1.K), b2.i, b2.j, b2.label)
            nb.extend([nb1, nb2, nb3, nb4])

        # Condition 3
        elif b1.i < b2.i and b1.j == b2.j:
            nb1 = Block(b1.K, b1.i, b2.i-1, b1.label[:b2.i-1-b1.i+1])
            nb2 = Block(b1.K, b2.i, b2.j, b2.label)
            nb3 = Block(set(b2.K)-set(b1.K),b2.i, b2.j, b2.label)
            nb.extend([nb1, nb2, nb3])

        # Condition4
        elif b1.i < b2.i and b2.i < b1.j and b1.j < b2.j:
            left_label   = b1.label[:b2.i-1-b1.i+1]
            center_label = b2.label[:b1.j-b2.i+1]
            right_label  = b2.label[b1.j-b2.i+1:]
            nb1 = Block(b1.K, b1.i, b2.i-1, left_label)
            nb2 = Block(set(b1.K)-set(b2.K), b2.i, b1.j, center_label)
            nb3 = Block(set(b1.K).intersection(b2.K), b2.i, b1.j, center_label)
            nb4 = Block(set(b2.K)-set(b1.K), b2.i, b1.j, center_label)
            nb5 = Block(b2.K, b1.j+1, b2.j, right_label)
            nb.extend([nb1, nb2, nb3, nb4, nb5])

    return nb