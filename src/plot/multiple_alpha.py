from pathlib import Path
from typing import List, Union, Optional
from .vertical_blocks import vertical_blocks

_Path = Union[str, Path]

def multiple_alpha(
                alphas: List[int],   
                path_msa: _Path, 
                path_blocks: _Path, 
                path_submsas: Optional[_Path] = None, 
                path_img: Optional[_Path] = None, 
                alpha: int=1, 
                rows_plot: int=100,
                ):
    
    imgs = []
    for alpha in alphas:
        imgs.append(
            vertical_blocks(
                path_msa, 
                path_blocks, 
                path_submsas, 
                path_img, 
                alpha, 
                rows_plot,
            )
        )    

    