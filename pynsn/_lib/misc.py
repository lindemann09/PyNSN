"""
Draw a random number from a beta dirstribution
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import sys
from collections import OrderedDict
from typing import Any

import numpy as np


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


def dict_to_csv(dictionary, variable_names=False, dict_of_lists=False):
    d = OrderedDict(dictionary.items())
    rtn = ""
    if variable_names:
        rtn += ",".join(d.keys()) + "\n"

    if dict_of_lists:
        prop_np = np.asarray(list(d.values())).T  # list is requires in PY3
        for x in prop_np:
            rtn += ", ".join(map(str, x)) + "\n"
    else:
        rtn += ",".join(map(str, d.values())) + "\n"

    return rtn


def key_value_format(key: str, value: Any) -> str:
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
        __IPYTHON__
        return True
    except NameError:
        pass

    is_idle = "idlelib.run" in sys.modules
    # ps2 is only defined in interactive mode
    return is_idle or hasattr(sys, "ps2")
