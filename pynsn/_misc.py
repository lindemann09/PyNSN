"""
Draw a random number from a beta distribution
"""

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import sys
from collections import OrderedDict
from typing import Any, List, Sequence, Union

import numpy as np
import orjson
from numpy.typing import ArrayLike, NDArray

IntOVector = Union[int, Sequence[int], NDArray[np.integer]]


def formated_json(d: dict, indent: int = 2) -> str:
    json_str = orjson.dumps(
        d, option=orjson.OPT_SERIALIZE_NUMPY).decode("utf-8")
    if indent < 1:
        return json_str
    rtn = ""
    i = 0
    block_newline = False
    for x in json_str:
        if x == "{":
            i = i + indent
            x = x + "\n" + " " * i
        elif x == "," and not block_newline:
            x = x + "\n" + " " * (i-1)  # -1, because a space follows
        elif x == "}":
            i -= indent
            x = "\n" + " " * i + x
        elif x == "[":
            block_newline = True
        elif x == "]":
            block_newline = False

        rtn += x

    return rtn


def delete_elements(lst: List, index: IntOVector) -> List:
    """delete multiple elements from list"""
    rtn = []

    if isinstance(index, (np.integer, int)):
        index = (index,)  # make iterable
    for c, s in enumerate(lst):
        if c not in index:
            rtn.append(s)
    return rtn


def join_dict_list(list_of_dicts):
    """make a dictionary of lists from a list of dictionaries"""
    rtn = OrderedDict()
    for d in list_of_dicts:
        for k, v in d.items():
            if k in rtn:
                rtn[k].append(v)
            else:
                rtn[k] = [v]
    return rtn


def key_value_format(key: str, value: Any) -> str:
    if isinstance(value, int):
        v = f"{value:14}"  # try rounding
    else:
        try:
            v = f"{value:14.2f}"  # try rounding
        except (ValueError, TypeError):
            v = f"{str(value):>14}"

    return f"{key:<24}{v}"


def dict_to_text(the_dict: dict) -> str:
    rtn = ""
    for k, v in the_dict.items():
        if len(rtn) == 0:
            rtn += "-"
        else:
            rtn += " "
        rtn += key_value_format(k, v) + "\n"
    return rtn.rstrip()


def is_interactive_mode():
    """Returns if Python is running in interactive mode (such as IDLE or IPthon)

    Returns
    -------
    interactive_mode : boolean

    """
    try:
        __IPYTHON__  # type: ignore
        return True
    except NameError:
        pass

    is_idle = "idlelib.run" in sys.modules
    # ps2 is only defined in interactive mode
    return is_idle or hasattr(sys, "ps2")


def polar2cartesian(polar: ArrayLike) -> NDArray:
    """polar has to be an 2D-array representing polar coordinates (radius, angle)"""
    polar = np.atleast_2d(polar)
    return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                    polar[:, 0] * np.sin(polar[:, 1])]).T


def cartesian2polar(xy: ArrayLike,
                    radii_only: bool = False) -> NDArray[np.floating]:
    """polar coordinates (radius, angle)

    if only radii required you may consider radii_only=True for faster
    processing

    xy has to be a 2D array
    """
    xy = np.atleast_2d(xy)
    rtn = np.hypot(xy[:, 0], xy[:, 1])
    if not radii_only:
        # add angle column
        rtn = np.array([rtn, np.arctan2(xy[:, 1], xy[:, 0])]).T
    return rtn
