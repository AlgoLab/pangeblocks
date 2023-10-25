"""Not needed anymore"""
from src.blocks import Block

def block_between_vertical_blocks(horizontal_block: Block, vertical_block1: Block, vertical_block2: Block):
    """If the horizontal block spans both vertical blocks, 
    return the block in the middle"""
    vb1, vb2 = list(sorted([vertical_block1, vertical_block2], key=lambda block: block.start))
    hb = horizontal_block

    if (hb.start <= vb1.start and vb1.end <= hb.end) and \
        (hb.start <= vb2.start and vb2.end <= hb.end) and \
        set(hb.K).issubset(vb1.K) and set(hb.K).issubset(vb1.K):
        
        start = vb1.end+1
        end   = vb2.start-1
        label = hb.label[start-hb.start : end-hb.start + 1]
        return Block(K=hb.K, start=start, end=end, label=label)

# Example
# >>> hb = Block([1,2,3],0,9,"ACGTACGTAC")
# >>> vb1 = Block([0,1,2,3,4,5],2,3,"GT")
# >>> vb2 = Block([0,1,2,3,4,5],7,8,"TA")
# >>> block_between_vertical_blocks(hb,vb1,vb2)

from collections import defaultdict
import logging
logging.basicConfig(level=logging.INFO,
                    format='[Block between vertical blocks]%(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')

def decomposition(vertical_blocks: list(Block), n_seqs: int, n_cols: int, msa):
    """Identify maximal blocks spanning two vertical blocks,
    and return a list with decomposed blocks, those that are in between two vertical
    blocks

    Args:
        input_blocks (list): _description_
        n_seqs (int): _description_
    """    
    input_blocks=vertical_blocks
    n_blocks = len(input_blocks)
    
    # covering_by_position is a dictionary with key (r,c) and value the list
    # of indices of the blocks that include the position (r,c)
    #
    # to speed up the process, we iterate over the blocks and append the
    # block to all the positions it covers
    covering_by_position = defaultdict(list)
    left_maximal_vertical_blocks = {}

    # We keep only maximal vertical blocks and we achieve that in two
    # phases:
    # 1. for each beginning column, we keep the block with the largest
    #    end column
    # 2. we do a second scan of vertical blocks and, for each end
    #    column, we keep the block with the smallest beginning column
    logging.info("loop: left maximal vertical blocks")
    for idx, block in enumerate(input_blocks):
        # logging.debug("input block: %s", block.str())
        if len(block.K) == n_seqs:
            if not block.start in left_maximal_vertical_blocks or block.end > left_maximal_vertical_blocks[block.start]["end"]:
                left_maximal_vertical_blocks[block.start] = {
                    "idx": idx, "end": block.end}

    logging.info("loop: vertical blocks")
    vertical_blocks = {}
    for begin, block in left_maximal_vertical_blocks.items():
        if not block["end"] in vertical_blocks or begin < vertical_blocks[block["end"]]["begin"]:
            vertical_blocks[block["end"]] = {
                "idx": block["idx"], "begin": begin}
    logging.info("vertical blocks: %s", vertical_blocks)

    # We keep a set of all columns that are covered by a vertical block:
    # those columns will not be involved in the U[r,c] variables and in the
    # corresponding covering constraints
    logging.info("loop: covered by vertical blocks")
    covered_by_vertical_block = set()
    for begin, item in vertical_blocks.items():
        block = input_blocks[item['idx']]
        # logging.info("Vertical block: %s, block.str()")
        for col in range(block.start, block.end + 1):
            covered_by_vertical_block.add(col)
    logging.info(
        "Covered by vertical blocks: %s", covered_by_vertical_block)
    logging.info(
        "No. covered by vertical blocks: %s out of %s" % (len(covered_by_vertical_block), self.n_cols))


    # msa_positions is a list of all positions (r,c) that are required to be
    # covered. We exclude the positions covered by vertical blocks, since
    # they will be guaranteed to be covered, as an effect of the fact that
    # the corresponding C variables will be set to 1
    logging.debug("generating msa positions")
    msa_positions = [(r, c) for r in range(n_seqs)
                        for c in set(range(n_cols)) - covered_by_vertical_block]

    # We compute a dictionary called zones with key the column and value the
    # zone_id, that is a progressive id of the region, where each region is
    # a set of consecutive columns that are either disjoint or included in a vertical block.
    logging.info("zones")
    current_zone = 0
    zone = {0: 0}
    zone_boundaries = {0: {'start': 0, 'end': -1}}
    # We keep also a dictionary called zone_boundaries with key the zone_id
    # and values the first and last column of the zone
    for col in range(1, n_cols):
        if (col in covered_by_vertical_block) != (col - 1 in covered_by_vertical_block):
            # logging.debug("Boundary at %s: %s-%s", col, (col in covered_by_vertical_block),
            #               (col - 1 in covered_by_vertical_block))

            zone_boundaries[current_zone]['end'] = col - 1
            current_zone += 1
            zone_boundaries[current_zone] = {'start': col, 'end': -1}
        zone[col] = current_zone
    zone_boundaries[current_zone]['end'] = n_cols - 1


    # We compute the set unsplit_blocks, corresponding to the
    # blocks that are disjoint from or contained in vertical blocks.
    # The blocks that we will encode with a C variable in the
    # ILP correspond to the blocks that are disjoint from vertical blocks or
    # are vertical blocks themselves.
    logging.debug("zones\n%s", zone)
    logging.debug("boundaries:\n%s", zone_boundaries)
    # We start with all vertical blocks, since they will be encoded with a C variable
    private_blocks = {
        (input_blocks[vertical_blocks[k]['idx']].start,
            input_blocks[vertical_blocks[k]['idx']].end,
            frozenset(input_blocks[vertical_blocks[k]['idx']].K)): input_blocks[vertical_blocks[k]['idx']] for k in vertical_blocks}
    logging.debug("private vertical blocks: %s", private_blocks)

    for idx, block in enumerate(input_blocks):
        logging.debug("Analyzing %s out of %s. Size=%sx%s. Block=%s %s %s",
                        idx+1, n_blocks, len(block.K), block.end - block.start + 1, block.start, block.end, block.K)
        # Compute the zone of the boundaries of the block
        (zone_start, zone_end) = (zone[block.start], zone[block.end])
        logging.debug("Zones: %s-%s" % (zone_start, zone_end))
        # current_columns is the set of columns covered by the current block
        if (zone_start == zone_end):
            # Since the current block is contained in a single zone, it is
            # included in a vertical block or disjoint from vertical blocks.
            # In both cases, we encode it with a C variable
            logging.debug("unsplit block: %s %s %s" %
                            (block.start, block.end, block.K))
            # If the block is included in a vertical block, we can
            # discard it.
            # Notice that all maximal vertical blocks have been added to
            # private blocks before this for loop.
            if block.start not in covered_by_vertical_block:
                private_blocks[(block.start, block.end,
                                frozenset(block.K))] = block
                logging.debug("Added %s -> %s", (block.start, block.end,
                                                frozenset(block.K)), (block.start, block.end, block.K))
            # logging.debug("disjoint vertical block: %s", block.str())
        else:
            # The current block overlaps the vertical blocks.
            # We need to decompose it into parts that are disjoint from the
            # vertical blocks.
            logging.info("splitting block: %s %s %s" %
                            (block.start, block.end, block.K))

            # logging.debug("block %s is not disjoint", block.str())
            # logging.debug("zone_start: %s, zone_end: %s",
            #               zone_start, zone_end)
            # logging.debug("covered_by_vertical_block: %s",
            #               covered_by_vertical_block)

            # logging.debug(
            #     "block %s intersects with at least two vertical blocks", block.str())

            for zone_id in range(zone_start, zone_end + 1):
                begin, end = max(zone_boundaries[zone_id]['start'],block.start), min(zone_boundaries[zone_id]['end'], block.end)
                if begin not in covered_by_vertical_block:
                    # the current zone is not covered by a vertical block
                    if block.end < end:
                        # this is the last zone of the block, so the
                        # last position is the end of the block, not the
                        # end of the zone
                        end = block.end
                    label = str(msa[block.K[0], begin:end+1].seq)
                    new_block = Block(block.K, begin, end, label)
                    # check correct label of the block for all rows
                    for r in new_block.K:
                        for pos_block, c in enumerate(range(new_block.start, new_block.end + 1)):

                            # sanity check of labels
                            r = int(r)
                            char_msa = msa[r,c].upper()
                            char_block = new_block.label[pos_block].upper() 
                            if char_msa != char_block:
                                logging.info("incorrect new_block covering position (%s,%s): %s" % (r,c,new_block.str()))
                                logging.info("label of new_block should be %s", str(msa[r, new_block.start:new_block.end + 1].seq))


                    logging.debug("considering adding: %s %s %s" %
                                    (new_block.start, new_block.end, new_block.K))
                    private_blocks[(new_block.start, new_block.end,
                                    frozenset(new_block.K))] = new_block
                    logging.debug("Added %s -> %s", (new_block.start, new_block.end,
                                                    frozenset(new_block.K)),  (new_block.start, new_block.end, new_block.K))

                    # logging.debug("Adding private block: %s to %s" % (new_block.str(), private_blocks))
    logging.debug("Private blocks: %s", private_blocks)

    # Remove duplicates
    good_blocks = list(private_blocks.values())
    vertical_blocks = [idx for idx, block in enumerate(good_blocks) if (
        zone[block.start] == zone[block.end]) and block.start in covered_by_vertical_block]
    del private_blocks
    logging.info(
        f"Total number of distinct blocks: {len(good_blocks)}")
    logging.info("Vertical blocks: %s", len(vertical_blocks))
    
    return good_blocks
