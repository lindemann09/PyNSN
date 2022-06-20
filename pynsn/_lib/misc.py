"""
Draw a random number from a beta dirstribution
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from collections import OrderedDict
import numpy as np

try:
    from math import log2
except:
    from math import log

    log2 = lambda x: log(x, 2)


def is_base_string(s):
    return isinstance(s, (str, bytes))

def is_unicode_string(s):
    return isinstance(s, str)

def is_byte_string(s):
    return isinstance(s, bytes)

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
        feat_np = np.array(list(d.values())).T  # list is requires in PY3
        for x in feat_np:
            rtn += ", ".join(map(lambda s: str(s), x)) + "\n"
    else:
        rtn += ",".join(map(lambda s: str(s), d.values())) + "\n"

    return rtn

def numpy_vector(x):
    """helper function:
    make an numpy vector from any element (list, arrays, and single data (str, numeric))
    """

    x = np.array(x)
    if x.ndim == 1:
        return x
    elif x.ndim == 0:
        return x.reshape(1)  # if one element only, make a array with one element
    else:
        return x.flatten()

def numpy_array_2d(two_d_data):
    """ensures well shaped to 2d numpy array"""
    rtn = np.array(two_d_data)
    if rtn.ndim == 1 and len(rtn) == 2:
        rtn = rtn.reshape((1, 2))
    if rtn.ndim != 2:
        raise ValueError("Bad shaped data: xy must be pair of xy-values or a list of xy-values")
    return rtn

def numpy_round2(array, decimals, int_type=np.int32):
    """rounds and changes to int type if decimals == 0"""
    array = np.round(array, decimals=decimals)
    if decimals == 0:
        return array.astype(int_type)
    else:
        return array


def is_all_larger(vector, standard=0):
    return sum(map(lambda x: x > standard, vector))==len(vector)

def is_all_smaller(vector, standard=0):
    return sum(map(lambda x: x < standard, vector))==len(vector)

def is_all_equal(vector):
    # returns true if all elements are equal
    return sum(map(lambda x: x==vector[0], vector))==len(vector)
