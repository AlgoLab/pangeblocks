from pydantic import root_validator, validator
from pydantic.dataclasses import dataclass as pydataclass
from dataclasses import dataclass
from .positional_string import PositionalString
from typing import Optional

@pydataclass(eq=True, frozen=True)
class Block:
    "class for keeping track a block"
    K: tuple
    start: int
    end: int
    label: str

    def to_positional_string(self,) -> PositionalString:
        return PositionalString(self.label, self.start, self.end)

    def str(self):
        return "%s,%s,%s,%s" % (self.K,self.start,self.end,self.label)

    def len(self):
        "Length of the blocks, number of columns it spans"
        return self.end-self.start+1
    
    @root_validator
    def check_len_label(cls, values):
        "length of string must be equal to end-start+1"
        label = values.get("label")
        start = values.get("start")
        end = values.get("end")

        if len(label) != end-start+1:
            raise ValueError("len of the string 'label' must match the distance between columns 'start' and 'end'")
        return values

    @root_validator
    def check_start_less_or_equal_end(cls, values):
        "start <=end"
        start = values.get("start")
        end = values.get("end")

        if start > end:
            raise ValueError(f"start must be <= than end, start={start} and end={end}")
        return values

    @validator("K")
    def sort_K(cls, v):
        "sort values in K"
        return tuple(sorted(v))

@pydataclass(eq=True, frozen=True)
class LightBlock:
    "class for keeping track a block without the label"
    K: tuple
    start: int
    end: int
        
    def str(self):
        return "%s,%s,%s" % (self.K,self.start,self.end)

    def len(self):
        return self.end-self.start+1

    def nrows(self):
        return len(self.K)
    
    def ncols(self):
        return self.end-self.start+1
    
    def ncells(self):
        return self.nrows()*self.ncols()
    
    @root_validator
    def check_start_less_or_equal_end(cls, values):
        "start <=end"
        start = values.get("start")
        end = values.get("end")

        if start > end:
            raise ValueError(f"start must be <= than end, start={start} and end={end}")
        return values

    @validator("K")
    def sort_K(cls, v):
        "sort values in K"
        return tuple(sorted(v))
    
    def __getitem__(self, item):
        if item==0:
            return self.K
        elif item==1 or item==-2: 
            return self.start 
        elif item==2 or item == -1: 
            return self.end

    # def __hash__(self) -> int:
    #     pass

    # def __eq__(self, __value: object) -> bool:
    #     pass