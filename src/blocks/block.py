from pydantic import root_validator, validator
from pydantic.dataclasses import dataclass
from ..positional_strings import PositionalString

@dataclass(eq=True, frozen=True)
class Block:
    "class for keeping track a block"
    K: tuple
    i: int
    j: int
    label: str

    def to_positional_string(self,) -> PositionalString:
        return PositionalString(self.label, self.i, self.j)

    def str(self):
        return f"{self.K},{self.i},{self.j},{self.label}"

    @root_validator
    def check_len_label(cls, values):
        "length of string must be equal to j-i+1"
        label = values.get("label")
        i = values.get("i")
        j = values.get("j")

        if len(label) != j-i+1:
            raise ValueError("len of the string 'label' must match the distance between columns 'i' and 'j'")
        return values

    @root_validator
    def check_i_less_or_equal_j(cls, values):
        "i <=j"
        i = values.get("i")
        j = values.get("j")

        if i > j:
            raise ValueError(f"i must be <= than j, i={i} and j={j}")
        return values

    @validator("K")
    def sort_K(cls, v):
        "sort values in K"
        return tuple(sorted(v))