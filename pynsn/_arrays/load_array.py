"""loading array
"""

import json
from typing import Union
from .dot_array import DotArray
from .rect_array import RectangleArray


def load_array(filename: str) -> Union[DotArray, RectangleArray, None]:
    """Loading json array file

    Args:
        filename: the file name

    Returns:
        Object Array
    """

    with open(filename, 'r') as fl:
        d = json.load(fl)

    arr_type = d["type"]
    if arr_type == "DotArray":
        return DotArray.from_dict(d)
    elif arr_type == "RectangleArray":
        return RectangleArray.from_dict(d)
    else:
        raise RuntimeError(f"Unknown array type {arr_type}")
