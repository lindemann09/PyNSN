"""
Draw a random number from a beta dirstribution
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from collections import OrderedDict
from tkinter import N
import numpy as np
import sys
from numpy.typing import NDArray, ArrayLike


def make_csv(xy, size_data_dict, attributes=None,
             array_hash=None, make_variable_names=True):
    """Not for User
    makes csv for Arrays with object of different size information
    size_data_dict: keys = variable names (e.g. width, height),
                    values vector of size data
    """
    rtn = ""
    if make_variable_names:
        if array_hash:
            rtn += "hash,"
        rtn += "x,y," + ",".join(size_data_dict.keys()) + ","
        if attributes is not None:
            rtn += "attribute,"
        rtn = rtn[:-1] + "\n"  # replace comma

    size_data = np.array(list(size_data_dict.values())).T
    if attributes is None:
        attribute_vector = [None] * len(xy)  # to have something to loop
    else:
        attribute_vector = attributes

    for pos, size, attr in zip(xy, size_data, attribute_vector):
        if array_hash:
            rtn += "{0},".format(array_hash)
        rtn += "{},{},".format(pos[0], pos[1])
        for s in size:
            rtn += "{},".format(s)
        if attributes is not None:
            rtn += "{},".format(attr)
        rtn = rtn[:-1] + "\n"  # replace comma
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


def dict_to_csv(dictionary, variable_names=False, dict_of_lists=False):
    d = OrderedDict(dictionary.items())
    rtn = ""
    if variable_names:
        rtn += ",".join(d.keys()) + "\n"

    if dict_of_lists:
        prop_np = np.asarray(list(d.values())).T  # list is requires in PY3
        for x in prop_np:
            rtn += ", ".join(map(lambda s: str(s), x)) + "\n"
    else:
        rtn += ",".join(map(lambda s: str(s), d.values())) + "\n"

    return rtn


def numpy_vector(x):
    """helper function:
    make an numpy vector from any element (list, _arrays, and single data (str, numeric))
    """

    x = np.array(x)
    if x.ndim == 1:
        return x
    elif x.ndim == 0:
        return x.reshape(1)  # if one element only, make a array with one element
    else:
        return x.flatten()


def numpy_array_2d(two_d_data: ArrayLike):
    """ensures well shaped to 2d numpy array"""
    rtn = np.asarray(two_d_data)
    if rtn.ndim == 1 and len(rtn) == 2:
        rtn = rtn.reshape((1, 2))
    if rtn.ndim != 2:
        raise ValueError("Bad shaped data: xy must be pair of xy-values or a list of xy-values")
    return rtn


def numpy_round2(array: NDArray[np.float_], decimals: int,
                 int_type: type =np.int32) -> NDArray:
    """rounds and changes to int type if decimals == 0"""
    array = np.round(array, decimals=decimals)
    if decimals == 0:
        return array.astype(int_type)
    else:
        return array


def is_all_equal(vector):
    # returns true if all elements are equal
    return len(np.unique(np.asarray(vector))) == 1


def dict_to_text(the_dict, col_a=22, col_b=14,
                 spacing_char=" "):

    rtn = ""
    for k, v in the_dict.items():
        if len(rtn) == 0:
            key_str = "- " + k
        else:
            key_str = "  " + k

        value = "{}\n".format(v)
        len_col_b = col_b - len(value)
        if len_col_b < 2:
            len_col_b = 2
        rtn += key_str + (spacing_char * (col_a - len(key_str))) + \
               (" " * len_col_b) + value

    return rtn.rstrip()


def is_interactive_mode():
    """Returns if Python is running in interactive mode (such as IDLE or IPthon)

    Returns
    -------
    interactive_mode : boolean

    """
    # ipython?
    try:
        __IPYTHON__
        return True
    except NameError:
        pass

    is_idle = "idlelib.run" in sys.modules
    # ps2 is only defined in interactive mode
    return is_idle or hasattr(sys, "ps2")


def triu_nan(m, k=0):
    """helper function
    upper triangular but nan instead of zeros (as in numpy's original function,
    see docu numpy.triu)
    """
    return m + np.tril(np.full(m.shape, np.nan), k=k - 1)
