from pydantic import (
    validate_arguments,
    # ValidationError,
    # validator,
)

from dataclasses import dataclass

@validate_arguments
@dataclass
class PositionalString:
    b: str
    i: int
    j: int

    def __len__(self,) -> int:
        return self.j-self.i + 1

    # TODO: validator to check len(b) == j-i+1