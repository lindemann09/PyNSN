from copy import copy as _copy
import itertools as _itertools

_DEFAULT_SPACING_PRECISION = 0.0001

# @total_ordering
class _BaseFeatureType(object):
    """"""
    label = "undefined"

    def __init__(self, value=None):
        self.value = value

    @property
    def type(self):
        return type(self).__name__

    def __add__(self, other):
        rtn = _copy(self)
        rtn.value += other.value
        return rtn

    def __sub__(self, other):
        rtn = _copy(self)
        rtn.value += other.value
        return rtn

    def __mul__(self, other):
        rtn = _copy(self)
        rtn.value *= other.value
        return rtn

    def __truediv__(self, other):
        rtn = _copy(self)
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


class _SizeRelatedFeatureType(_BaseFeatureType):

    @property
    def dependencies(self):
        return SIZE_FEATURES

class _SpaceRelatedFeatureType(_BaseFeatureType):

    def __init__(self, value=None, spacing_precision=_DEFAULT_SPACING_PRECISION):
        _BaseFeatureType.__init__(self, value)
        self.spacing_precision = spacing_precision

    @property
    def dependencies(self):
        return SPACE_FEATURES

    def as_dict(self):
        d = _BaseFeatureType.as_dict(self)
        d["spacing_precision"] = self.spacing_precision
        return (d)


class LogSize(_SizeRelatedFeatureType):
    label = "Log Size"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.logSize


class LogSpacing(_SpaceRelatedFeatureType):
    label = "Log Spacing"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.logSpacing


class TotalSurfaceArea(_SizeRelatedFeatureType):
    """"""
    label = "Total surface area"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.total_surface_area


class ItemDiameter(_SizeRelatedFeatureType):
    """"""
    label = "Mean item diameter"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.mean_item_diameter

class ItemSurfaceArea(_SizeRelatedFeatureType):
    """"""
    label = "Mean item surface area"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.mean_item_surface_area


class ItemPerimeter(_SizeRelatedFeatureType):
    """"""
    label = "Mean item perimeter"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.mean_item_perimeter

class TotalPerimeter(_SizeRelatedFeatureType):
    """"""
    label = "Total perimeter"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.total_perimeter


class Sparsity(_SpaceRelatedFeatureType):
    """"""
    label = "Sparsity"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.sparsity


class FieldArea(_SpaceRelatedFeatureType):
    """"""
    label = "Field area"

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.field_area


class Coverage(_BaseFeatureType):
    """"""
    label = "Coverage"

    def __init__(self, value=None, match_ratio_fieldarea2totalarea=0.5, spacing_precision=None):

        _BaseFeatureType.__init__(self, value)
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

    def adapt_value(self, reference_dot_array):
        self.value = reference_dot_array.feature.converage

    @property
    def dependencies(self):
        dep = [Coverage]
        if self.match_ratio_fieldarea2totalarea < 1:  # area involved
            dep.extend(TotalSurfaceArea().dependencies)
        if self.match_ratio_fieldarea2totalarea > 0:  # convex hull involved
            dep.extend(FieldArea().dependencies)
        return dep

    def as_dict(self):
        d = _BaseFeatureType.as_dict(self)
        d["match_ratio_convhull2area"] = self.match_ratio_fieldarea2totalarea
        d["spacing_precision"] = self.spacing_precision
        return (d)


ALL_VISUAL_FEATURES = (ItemDiameter, TotalSurfaceArea, TotalPerimeter,
                       FieldArea, Coverage, LogSpacing, LogSize, Sparsity,
                       ItemSurfaceArea, ItemPerimeter)

SIZE_FEATURES = (ItemDiameter, ItemSurfaceArea, ItemPerimeter,
                 TotalPerimeter, TotalSurfaceArea, LogSize)

SPACE_FEATURES = (LogSpacing, FieldArea, Sparsity) #TODO wha about Coverage?



def check_feature_list(feature_list, check_set_value=False):
    """helper function
    raises TypeError or Runtime errors if checks fail
    * type check
    * dependency check
    * is value defined
    """
    for x in feature_list:
        if not isinstance(x, _BaseFeatureType):
            raise TypeError("Parameter is not a continuous properties or a " + \
                            "list of continuous properties") #FIXME labels
            # continious property or visual feature
        elif check_set_value and x.value is None:
            raise RuntimeError("Value of continuous property {} is not defined.".format(
                type(x).__name__))

    # check dependencies
    for a, b in _itertools.combinations(feature_list, 2):
        if a.is_dependent(b):
            raise RuntimeError("Incompatible properties to match: {} & {}".format(
                type(a).__name__, type(b).__name__))