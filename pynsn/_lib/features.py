from copy import copy
import itertools

DEFAULT_SPACING_PRECISION = 0.0001

# @total_ordering
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
    def dependencies(self):
        return ()

    def is_dependent(self, other):
        return type(other) in self.dependencies or \
               type(self) in other.dependencies

    def as_dict(self):
        return {"type": self.type,
                "value": self.value}


class _SizeRelatedFeature(_BaseFeature):

    @property
    def dependencies(self):
        return (ItemDiameter, ItemSurfaceArea, ItemPerimeter, TotalPerimeter, TotalSurfaceArea, LogSize)


class _SpaceRelatedFeature(_BaseFeature):

    def __init__(self, value=None, spacing_precision=DEFAULT_SPACING_PRECISION):
        _BaseFeature.__init__(self, value)
        self.spacing_precision = spacing_precision

    @property
    def dependencies(self):
        return (LogSpacing, FieldArea, Sparsity)

    def as_dict(self):
        d = _BaseFeature.as_dict(self)
        d["spacing_precision"] = self.spacing_precision
        return (d)


class LogSize(_SizeRelatedFeature):
    long_label = "Log Size"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_logSize


class LogSpacing(_SpaceRelatedFeature):
    long_label = "Log Spacing"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_logSpacing


class TotalSurfaceArea(_SizeRelatedFeature):
    """"""
    long_label = "Total surface area"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_total_surface_area


class ItemDiameter(_SizeRelatedFeature):
    """"""
    long_label = "Mean item diameter"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_mean_item_diameter

class ItemSurfaceArea(_SizeRelatedFeature):
    """"""
    long_label = "Mean item surface area"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_mean_item_surface_area


class ItemPerimeter(_SizeRelatedFeature):
    """"""
    long_label = "Mean item perimeter"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_mean_item_perimeter

class TotalPerimeter(_SizeRelatedFeature):
    """"""
    long_label = "Total perimeter"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_total_perimeter


class Sparsity(_SpaceRelatedFeature):
    """"""
    long_label = "Sparsity"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_sparsity


class FieldArea(_SpaceRelatedFeature):
    """"""
    long_label = "Field area"

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.feature_field_area


class Coverage(_BaseFeature):
    """"""
    long_label = "Coverage"

    def __init__(self, value=None, match_ratio_fieldarea2totalarea=0.5, spacing_precision=None):

        _BaseFeature.__init__(self, value)
        self.match_ratio_fieldarea2totalarea = match_ratio_fieldarea2totalarea
        if spacing_precision is None:
            self.spacing_precision = FieldArea().spacing_precision
        else:
            self.spacing_precision = spacing_precision

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
    def dependencies(self):
        dep = [Coverage]
        if self.match_ratio_fieldarea2totalarea < 1:  # area involved
            dep.extend(TotalSurfaceArea().dependencies)
        if self.match_ratio_fieldarea2totalarea > 0:  # convex hull involved
            dep.extend(FieldArea().dependencies)
        return dep

    def as_dict(self):
        d = _BaseFeature.as_dict(self)
        d["match_ratio_convhull2area"] = self.match_ratio_fieldarea2totalarea
        d["spacing_precision"] = self.spacing_precision
        return (d)


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
