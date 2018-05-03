from __future__ import division

from copy import copy
import itertools

#@total_ordering
class _BaseFeature(object):
    """"""
    long_label = "undefined"

    def __init__(self, value=None):
        self.value = value

    @property
    def type(self):
        return type(self).__name__

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
        return {"type": self.type,
                "value": self.value}


class TotalSurfaceArea(_BaseFeature):
    """"""
    long_label = "Total surface area"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_total_surface_area

    @property
    def _dependencies(self):
        return [ItemDiameter, TotalPerimeter, TotalSurfaceArea]



class ItemDiameter(_BaseFeature):
    """"""
    long_label = "Mean item diameter"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_item_diameter

    @property
    def _dependencies(self):
        return [ItemDiameter, TotalPerimeter, TotalSurfaceArea]


class TotalPerimeter(_BaseFeature):
    """"""
    long_label = "Total perimeter"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_total_perimeter

    @property
    def _dependencies(self):
        return [ItemDiameter, TotalPerimeter, TotalSurfaceArea]


class FieldArea(_BaseFeature):
    """"""
    long_label = "Field area"

    def __init__(self, value=None, match_presision=0.01):
        _BaseFeature.__init__(self, value)
        self.match_presision = match_presision

    @property
    def _dependencies(self):
        return [FieldArea]

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_field_area

    def as_dict(self):
        d = _BaseFeature.as_dict(self)
        d["match_presision"] = self.match_presision
        return(d)


class Coverage(_BaseFeature):
    """"""
    long_label = "Coverage"

    def __init__(self, value=None, match_ratio_fieldarea2totalarea=0.5, convex_hull_precision=None):

        _BaseFeature.__init__(self, value)
        self.match_ratio_fieldarea2totalarea = match_ratio_fieldarea2totalarea
        if convex_hull_precision is None:
            self.convex_hull_precision = FieldArea().match_presision
        else:
            self.convex_hull_precision = convex_hull_precision

    @property
    def match_ratio_fieldarea2totalarea(self):
        return self._match_ratio_fieldarea2totalarea

    @match_ratio_fieldarea2totalarea.setter
    def match_ratio_fieldarea2totalarea(self, v):
        self._match_ratio_fieldarea2totalarea = v
        if v < 0 or v > 1:
            raise RuntimeError("Match_ratio_convhull2area has to be between 0 and 1.")

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_converage

    @property
    def _dependencies(self):
        dep = [Coverage]
        if self.match_ratio_fieldarea2totalarea < 1:  # area involved
            dep.extend(TotalSurfaceArea()._dependencies)
        if self.match_ratio_fieldarea2totalarea > 0:  # convex hull involved
            dep.extend(FieldArea()._dependencies)
        return dep

    def as_dict(self):
        d = _BaseFeature.as_dict(self)
        d["match_ratio_convhull2area"] = self.match_ratio_fieldarea2totalarea
        d["convex_hull_precision"] = self.convex_hull_precision
        return(d)

## helper function
def check_feature_list(feature_list, check_set_value=False):
    """helper function
    raises TypeError or Runtime errors if checks fail
    * type check
    * dependency check
    * checks if value is defined
    """
    for x in feature_list:
        if not isinstance(x, _BaseFeature):
            raise TypeError("Parameter is not a continuous properties or a " + \
                            "list of continuous properties")
        elif check_set_value and x.value is None:
            raise RuntimeError("Value of continuous property {} is not defined.".format(
                type(x).__name__))

    # check dependencies
    for a, b in itertools.combinations(feature_list, 2):
        if a.is_dependent(b):
            raise RuntimeError("Incompatible properties to match: {} & {}".format(
                type(a).__name__, type(b).__name__))

def _dict_to_feature(d):
    d = copy(d)
    t = d["type"]
    del d["type"]
    if t == ItemDiameter.__name__:
        return ItemDiameter(**d)
    elif t == TotalSurfaceArea.__name__:
        return TotalSurfaceArea(**d)
    elif t == TotalPerimeter.__name__:
        return TotalPerimeter(**d)
    elif t == Coverage.__name__:
        print("l")
        return Coverage(**d)
    elif t == FieldArea.__name__:
        return FieldArea(**d)
    return None

