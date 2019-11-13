"""
Draw a random number from a beta dirstribution
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from collections import OrderedDict
import random
import numpy as np
from scipy import spatial

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

def polar2cartesian(polar):
    """polar is an 2d-array representing polar coordinates (radius, angle)"""
    polar = np.array(polar)
    return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                     polar[:, 0] * np.sin(polar[:, 1])]).T

def cartesian2polar(xy, radii_only=False):
    """polar coordinates (radius, angle)

    if only radii required you may consider radii_only=True for faster
    processing
    """
    xy = np.array(xy)
    radii = np.hypot(xy[:, 0], xy[:, 1])
    if radii_only:
        return radii
    else:
        return np.array([radii, np.arctan2(xy[:, 1], xy[:, 0])]).T


class EfficientConvexHull(object):
    """helper class for efficent (and not repeated) calulations convex hull
    convex hull not be recalulated
    """

    def __init__(self, xy):
        self._xy = xy
        self._scipy_ch_object = None

    @property
    def scipy_convex_hull(self):
        if self._scipy_ch_object is None:
            self._scipy_ch_object = spatial.ConvexHull(self._xy)

        return self._scipy_ch_object


    @property
    def indices(self):
        return self.scipy_convex_hull.vertices

    @property
    def xy(self):
        return self._xy[self.indices, :]


class EfficientConvexHullDots(EfficientConvexHull):
    """helper class for efficent (and not repeated) calulations convex hull"""

    def __init__(self, xy, diameters):
        EfficientConvexHull.__init__(self, xy=xy)
        self._diameters = diameters
        self._full_xy = None
        self._full_area = None

    @property
    def full_xy(self):
        """this convex_hull takes into account the dot diameter"""
        if self._full_xy is None:
            idx = self.scipy_convex_hull.vertices

            minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)  # center outer positions

            polar_centered = cartesian2polar(self._xy[idx, :] - center)
            polar_centered[:, 0] = polar_centered[:, 0] + (self._diameters[idx] / 2)
            self._full_xy = polar2cartesian(polar_centered) + center

        return self._full_xy

    @property
    def full_field_area(self):
        if self._full_area is None:
            self._full_area = spatial.ConvexHull(self.full_xy).volume

        return self._full_area
