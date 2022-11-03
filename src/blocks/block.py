from pydantic import validate_arguments
from dataclasses import dataclass
from ..positional_strings import PositionalString

@validate_arguments
@dataclass(eq=True,frozen=True)
class Block:
    "class for keeping track a block"
    K: tuple
    i: int
    j: int
    label: str

    def to_positional_string(self,) -> PositionalString: 
        return PositionalString(self.label, self.i, self.j)

    # TODO: validator to check len(b) == j-i+1