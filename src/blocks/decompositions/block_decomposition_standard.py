from .. import Block


import logging
logging.basicConfig(level=logging.INFO,
                    format='[Block Decomposition Row Maximal] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

def block_decomposition(block1: Block, block2: Block):
    """Decomposition of maximal blocks
    two cases according to lemma 1 in the manuscript
    """
    logging.debug("standard decomposition")
    # sort blocks by starting positions
    B = [block1, block2]
    l2, l1 = list(sorted(B, key=lambda l: (l.start,len(l.K))))

    # not empty intersection
    common_rows = list(set(l1.K).intersection(set(l2.K)))
    common_cols = list(set(range(l1.start,l1.end+1)).intersection(set(range(l2.start,l2.end+1))))
    
    if common_rows and common_cols: 

        logging.debug(f"l1: {l1}")
        logging.debug(f"l2: {l2}")

        # Lemma 1: [b1,e1] subset of [b2,e2] <=> K2 subset of K1 (l2 is the most one to the left)
        blocks = []
        if l2.start <= l1.start and l1.end <= l2.end: 

            if l2.start <= l1.start -1: blocks.append( Block(l2.K , l2.start, l1.start - 1))
            if l1.end + 1 <= l2.end:  blocks.append( Block(l2.K, l1.end + 1, l2.end))
        
            blocks.extend(
                [
                Block(set(l1.K).intersection(l2.K), l1.start, l1.end),
                Block(set(l1.K).difference(l2.K), l1.start, l1.end)
                ]
                )
        else:
            l1, l2 = list(sorted(B, key=lambda l: (l.start, len(l.K)))) # l1 es el de mas a la izquierda
            blocks = [
                Block(l1.K, l1.start, l2.start-1),
                Block(set(l1.K).difference(l2.K), l2.start, l1.end),
                Block(set(l1.K).intersection(l2.K), l2.start, l1.end),
                Block(set(l2.K).difference(l1.K), l2.start, l1.end),
                Block(l2.K, l1.end + 1, l2.end)
            ]
        
        return blocks