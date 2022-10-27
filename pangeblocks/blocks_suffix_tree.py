"""
Build blocks (K,i,j) for a MSA of length n
K: set of rows
i: start position (from [0,n-1])
j: end position (from [1,n])
"""
import numpy as np
from collections import defaultdict, namedtuple
from suffix_tree import Tree
from typing import Tuple, List, Union

ALPHABET ={"A","C","G","T","-"}
FIRST_CHAR = chr(100)
## ___ ___ ___ 
Block = namedtuple("Block",["K","i","j","seq"])
_Block = Tuple[List[Union[int,str]], int, int] # Block (K,i,j)

## Functions
def save_as_dot(tree, path_save: str = "graph.dot"):
    " Save tree as .dot file to generate plot with grahviz"
    with open("graph.dot", "w") as fp:
        fp.write(tree.to_dot())

def decorate_seq(seq: str):
    global FIRST_CHAR
    "add characters between elements of the original sequence"
    dec_seq=[]
    for j,s in enumerate(seq):
        dec_seq.extend([chr(ord(FIRST_CHAR)+j),s])
    return "".join(dec_seq)

def inverse_decorate_seq(dec_seq: str):
    global ALPHABET
    "extract only characters in the orignal alphabet"
    return "".join([c for c in list(dec_seq) if c in ALPHABET])

def is_vertical_block(block: _Block, size_msa: int) -> bool:
    K,i,j = block
    if len(K)== size_msa and i<=j:
        return True
    return False

def maximal_blocks_from(i, seqs) -> List[_Block]:
    "Maximal blocks of sequences from position i"
    # define tree for sequences starting at position i
    tree = Tree({num: decorate_seq(seq[i:]) for num, seq in enumerate(seqs)})

    # Find Blocks starting at position i
    potential_maxblocks=defaultdict(int)
    for C, path in sorted(tree.maximal_repeats()):
        seq, start, end= "".join(list(path.S[:-1])) , path.start, path.end
        seq_block = seq[start:end]

        # consider only blocks starting at position i (len>2 cause all start with FIRST_CHAR)
        if len(seq_block)>1 and seq_block.startswith(FIRST_CHAR): 
            potential_maxblocks[seq[start:end]] +=1


        
    # filter blocks (K,i,j) for a MSA of length n
    max_blocks = []
    for seq_block in potential_maxblocks.keys():
        K = []
        seq_block_alphabet = inverse_decorate_seq(seq_block)

        for seq_id, seq in enumerate(seqs):
            if seq[i:].startswith(seq_block_alphabet):
                K.append(seq_id)
        
        # only starting at position i
        if len(K)>1:
            max_blocks.append(Block(K,i,i+len(seq_block_alphabet)-1,seq_block_alphabet))
    
    return max_blocks

def is_maximal_block(block, previous_blocks) -> bool:
    "Compare new block against maximal blocks"
    for pblock in previous_blocks:
        # if not a maximal block return False 
        if set(block.K) == set(pblock.K) and block.j==pblock.j and block.i > pblock.i:
            return False
    return True

def map_coverage(blocks, size_msa, n_seqs):
    "fill a matrix with 1 if there exist a maximal block covering a position"
    coverage = np.zeros((n_seqs, size_msa))
    for block in blocks: 
        for k in block.K:
            coverage[k,block.i:block.j+1] = 1
    return coverage

