"""
Draw a random number from a beta dirstribution
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from collections import OrderedDict
import random
import numpy as np
from ._visual_features import Features

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

def random_beta(size, number_range, mean, std):
    """Draw from beta distribution defined by the
    number_range [a,b] and mean and standard distribution

    Resulting distribution has the defined mean and std

    for calculated shape parameters [alpha, beta] see `shape_parameter_beta`

    Parameter:
    ----------
    number_range : tuple (numeric, numeric)
        the range of the distribution
    mean: numeric
    std: numeric TODO

    Note:
    -----
    Depending on the position of the mean in number range the
    distribution is left or right skewed.

    """


    if std is None or number_range is None or std == 0:
        return np.array([mean]*size)

    alpha, beta = shape_parameter_beta(number_range=number_range,
                                       mean=mean,
                                       std=std)

    # NOTE: do not use np.random.beta, because it produces identical numbers for
    # different threads:
    dist = np.array([random.betavariate(alpha=alpha, beta=beta) \
                     for _ in range(size)])
    dist = (dist - np.mean(dist)) / np.std(dist) # z values
    return dist*std + mean

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


def is_all_larger(vector, standard=0):
    return sum(map(lambda x: x > standard, vector))==len(vector)

def is_all_smaller(vector, standard=0):
    return sum(map(lambda x: x < standard, vector))==len(vector)

def is_all_equal(vector):
    # returns true if all elements are equal
    return sum(map(lambda x: x==vector[0], vector))==len(vector)



def check_feature_list(feature_list):
    """helper function
    raises TypeError or Runtime errors if checks fail
    * type check
    * dependency check
    """

    size_occured = ""
    space_occured = ""
    error = "Incompatible properties to match: {} & {}"

    if not isinstance(feature_list, (tuple, list)):
        feature_list = [feature_list]

    for x in feature_list:
        if x not in Features.ALL_FEATURES:
            raise TypeError("Parameter is not a continuous feature or a " + \
                            "list of continuous properties")
            # continious property or visual feature

        if x in Features.SIZE_FEATURES:
            if len(size_occured)>0:
                raise RuntimeError(error.format(x, size_occured))
            else:
                size_occured = x

        if x in Features.SPACE_FEATURES:
            if len(space_occured)>0:
                raise RuntimeError(error.format(x, space_occured))
            else:
                space_occured = x

