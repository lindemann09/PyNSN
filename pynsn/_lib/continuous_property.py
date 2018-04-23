from __future__ import division
from copy import copy
import itertools
from functools import total_ordering


#@total_ordering
class _ContinuousProperty(object):
    """"""
    long_label = "undefined"

    def __init__(self, value=None):
        self.value = value

    #def __lt__(self, other):
    #    return self.value < other.value

    #def __eq__(self, other):
    #    return self.value == other.value

    #def __ne__(self, other):
    #    return self.value != other.value

    def __add__(self, other):
        rtn = copy(self)
        rtn.value += other.value
        return rtn

    def __sub__(self, other):
        rtn = copy(self)
        rtn.value += other.value
        return rtn

    def __mul__(self, other):
        rtn = copy(self)
        rtn.value *= other.value
        return rtn

    def __truediv__(self, other):
        rtn = copy(self)
        rtn.value /= other.value
        return rtn

    def __int__(self):
        return int(self.value)

    def __float__(self):
        return float(self.value)

    @property
    def _dependencies(self):
        return []

    def is_dependent(self, other):
        return type(other) in self._dependencies or \
               type(self) in other._dependencies

    def as_dict(self):
        return {"type": type(self).__name__,
                "value": self.value}


class SurfaceArea(_ContinuousProperty):
    """"""
    long_label = "Total surface area"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_total_surface_area

    @property
    def _dependencies(self):
        return [DotDiameter, Circumference, SurfaceArea]



class DotDiameter(_ContinuousProperty):
    """"""
    long_label = "Mean dot diameter"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_mean_dot_diameter

    @property
    def _dependencies(self):
        return [DotDiameter, Circumference, SurfaceArea]


class Circumference(_ContinuousProperty):
    """"""
    long_label = "Total circumference"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_total_circumference

    @property
    def _dependencies(self):
        return [DotDiameter, Circumference, SurfaceArea]


class ConvexHull(_ContinuousProperty):
    """"""
    long_label = "Convex hull area"

    def __init__(self, value=None, match_presision=0.01):
        _ContinuousProperty.__init__(self, value)
        self.match_presision = match_presision

    @property
    def _dependencies(self):
        return [ConvexHull]

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_convex_hull_area

    def as_dict(self):
        d = _ContinuousProperty.as_dict(self)
        d["match_presision"] = self.match_presision
        return(d)


class Density(_ContinuousProperty):
    """"""
    long_label = "Density"

    def __init__(self, value=None, match_ratio_convhull2area=0.5, convex_hull_precision=None):

        _ContinuousProperty.__init__(self, value)
        self.match_ratio_convhull2area = match_ratio_convhull2area
        if convex_hull_precision is None:
            self.convex_hull_precision = ConvexHull().match_presision
        else:
            self.convex_hull_precision = convex_hull_precision

    @property
    def match_ratio_convhull2area(self):
        return self._match_ratio_convhull2area

    @match_ratio_convhull2area.setter
    def match_ratio_convhull2area(self, v):
        self._match_ratio_convhull2area = v
        if v < 0 or v > 1:
            raise RuntimeError("Match_ratio_convhull2area has to be between 0 and 1.")

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_density

    @property
    def _dependencies(self):
        dep = [Density]
        if self.match_ratio_convhull2area < 1:  # area involved
            dep.extend(SurfaceArea()._dependencies)
        if self.match_ratio_convhull2area > 0:  # convex hull involved
            dep.extend(ConvexHull()._dependencies)
        return dep

    def as_dict(self):
        d = _ContinuousProperty.as_dict(self)
        d["match_ratio_convhull2area"] = self.match_ratio_convhull2area
        d["convex_hull_precision"] = self.convex_hull_precision
        return(d)

## helper function
def check_list_continuous_properties(lcp, check_set_value=False):
    """helper function
    raises TypeError or Runtime errors if checks fail
    * type check
    * dependency check
    * checks if value is defined
    """
    for x in lcp:
        if not isinstance(x, _ContinuousProperty):
            raise TypeError("Parameter is not a continuous properties or a " + \
                            "list of continuous properties")
        elif check_set_value and x.value is None:
            raise RuntimeError("Value of continuous property {} is not defined.".format(
                type(x).__name__))

    # check dependencies
    for a, b in itertools.combinations(lcp, 2):
        if a.is_dependent(b):
            raise RuntimeError("Incompatible properties to match: {} & {}".format(
                type(a).__name__, type(b).__name__))

def dict_to_property(d):
    d = copy(d)
    t = d["type"]
    del d["type"]
    if t == DotDiameter.__name__:
        return DotDiameter(**d)
    elif t == SurfaceArea.__name__:
        return SurfaceArea(**d)
    elif t == Circumference.__name__:
        return Circumference(**d)
    elif t == Density.__name__:
        print("l")
        return Density(**d)
    elif t == ConvexHull.__name__:
        return ConvexHull(**d)
    return None

