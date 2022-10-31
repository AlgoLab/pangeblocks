from ..blocks.block import Block
from .positional_string import PositionalString

def positional_string_from_block(block: Block) -> PositionalString:
    return PositionalString(block.label, block.i, block.j)