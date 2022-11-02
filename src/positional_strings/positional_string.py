from pydantic import validate_arguments
from dataclasses import dataclass

@validate_arguments
@dataclass
class PositionalString:
    b: str
    i: int
    j: int