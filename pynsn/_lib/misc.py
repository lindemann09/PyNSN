"""
Draw a random number from a beta dirstribution
"""

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import sys
from collections import OrderedDict
from typing import Any


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
        __IPYTHON__
        return True
    except NameError:
        pass

    is_idle = "idlelib.run" in sys.modules
    # ps2 is only defined in interactive mode
    return is_idle or hasattr(sys, "ps2")
