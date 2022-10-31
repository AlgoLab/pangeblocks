from dataclasses import dataclass

@dataclass
class Block:
    "class for keeping track a block"
    K: set
    i: int
    j: int
    label: str
