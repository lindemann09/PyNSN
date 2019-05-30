"""
Draw a random number from a beta dirstribution
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import sys
import random
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

# randomizing

random.seed()


def random_beta(number_range, parameter):
    """Draw from beta distribution defined by the
    number_range [a,b] and the shape parameters [alpha, beta]

    Parameter:
    ----------
    number_range : tuple (numeric, numeric)
        the range of the distribution
    parameter: tuple
        shape parameter (alpha, beta) of the distribution
        see parameter function()

    Note:
    -----
    Depending on the position of the mean in number range the
    distribution is left or right skewed.

    """

    return number_range[0] + (number_range[1] - number_range[0]) \
           * random.betavariate(alpha=parameter[0], beta=parameter[1])


def shape_parameter_beta(number_range, mean, std):
    """Returns alpha (p) & beta (q) parameter for the beta distribution
    http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm

    Parameter
    ---------
    number_range : tuple (numeric, numeric)
        the range of the distribution
    mean : numeric
        the distribution mean
    std : numeric
         the distribution standard deviation

    Returns
    -------
    parameter: tuple
        shape parameter (alpha, beta) of the distribution

    """

    if mean <= number_range[0] or mean >= number_range[1] or \
            number_range[0] >= number_range[1]:
        raise RuntimeError("Mean has to be inside the defined number range")
    f = float(number_range[1] - number_range[0])
    mean = (mean - number_range[0]) / f
    std = (std) / f
    x = (mean * (1 - mean) / std ** 2 - 1)
    return (x * mean, (1 - mean) * x)


def join_dict_list(list_of_dicts):
    """make a dictionary of list from list of dictionaries"""
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
    nut None will not be procesed and returns None"""

    if x is None:
        return None

    x = np.array(x)
    if x.ndim == 1:
        return x
    elif x.ndim == 0:
        return x.reshape(1)  # if one element only, make a array with one element
    else:
        return x.flatten()


def is_all_larger(vector, standard=0):
    return sum(map(lambda x: x > standard, vector))==len(vector)

def is_all_smaller(vector, standard=0):
    return sum(map(lambda x: x < standard, vector))==len(vector)

def polar2cartesian(polar):
    """polar is an 2d-array representing polar coordinates (radius, angle)"""
    polar = np.array(polar)
    return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                     polar[:, 0] * np.sin(polar[:, 1])]).T

def cartesian2polar(xy, radii_only=False):
    xy = np.array(xy)
    radii = np.hypot(xy[:, 0], xy[:, 1])
    if radii_only:
        return radii
    else:
        return np.array([radii, np.arctan2(xy[:, 1], xy[:, 0])]).T
