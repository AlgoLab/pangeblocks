from pydantic import root_validator
from pydantic.dataclasses import dataclass

@dataclass
class PositionalString:
    b: str
    i: int
    j: int

    def __len__(self,) -> int:
        return self.j-self.i + 1

    @root_validator
    def check_len_b(cls, values):
        "lenght of string must be equal to j-i+1"
        b = values.get("b")
        i = values.get("i")
        j = values.get("j")

        if len(b) != j-i+1:
            raise ValueError("len of the string 'b' must match the distance between columns 'i' and 'j'")
        return values

    @root_validator
    def check_i_less_or_equal_j(cls, values):
        "i <=j"
        i = values.get("i")
        j = values.get("j")

        if i > j:
            raise ValueError(f"i must be <= than j, i={i} and j={j}")
        return values
