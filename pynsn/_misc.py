"""
Draw a random number from a beta distribution
"""

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import json
import sys
from collections import OrderedDict
from typing import Any, List, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import defaults

IntOVector = int | Sequence[int] | NDArray[np.integer]


def formated_json(d: dict, indent: int, decimals: None | int) -> str:
    if decimals is None:
        decimals = defaults.JSON_FLOAT_ROUNDING
    x = dict_convert_numpy(d, decimals)
    json_str = json.dumps(x, indent=None)
    if indent < 1:
        return json_str
    rtn = ""
    i = 0
    block_newline = False
    for x in json_str:
        if x == "{":
            i = i + indent
            x = x + "\n" + " " * (i - 1)
        elif x == "," and not block_newline:
            x = x + "\n" + " " * (i - 1)  # -1, because a space follows
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


def dict_of_arrays(array_of_dicts: list):
    rtn = {}
    # Loop over each key in the first dictionary
    for key in array_of_dicts[0].keys():
        # Extract each key's values across all dictionaries
        rtn[key] = [d[key] for d in array_of_dicts]
    return rtn


def dict_convert_numpy(d: dict, decimals: None | int = None) -> dict:
    """converts all numpy arrays in a dictionary to lists of Python scalars
    and, optionally, tries to round floating point data.

    """

    rtn = {}
    for k, v in d.items():
        if isinstance(v, np.ndarray):
            if np.issubdtype(v.dtype, np.integer):
                rtn[k] = v.tolist()
            elif np.issubdtype(v.dtype, np.floating):
                if isinstance(decimals, int):
                    v = np.round(v, decimals=decimals)
                rtn[k] = v.astype(float).tolist()
            else:
                rtn[k] = v.tolist()
                # raise TypeError(f"Can't handle dtype={v.dtype}, "
                #                "only integer or floating arrays")

        elif isinstance(v, dict):
            rtn[k] = dict_convert_numpy(v, decimals=decimals)
        elif isinstance(v, np.floating):
            if isinstance(decimals, int):
                v = np.round(v, decimals=decimals)
            rtn[k] = float(v)
        elif isinstance(v, np.integer):
            rtn[k] = int(v)
        else:
            rtn[k] = v

    return rtn


def is_interactive_mode():
    """Returns if Python is running in interactive mode (such as IDLE or IPython)

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
    return np.array(
        [polar[:, 0] * np.cos(polar[:, 1]), polar[:, 0] * np.sin(polar[:, 1])]
    ).T


def cartesian2polar(xy: ArrayLike, radii_only: bool = False) -> NDArray[np.floating]:
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
