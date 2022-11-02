from pydantic import validate_arguments
from dataclasses import dataclass

@validate_arguments
@dataclass(eq=True,frozen=True)
class Block:
    "class for keeping track a block"
    K: tuple
    i: int
    j: int
    label: str
