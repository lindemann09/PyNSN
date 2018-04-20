from functools import total_ordering
import itertools

@total_ordering
class _ContinuousProperty(object):
    """"""
    _AREA = "TA"
    _DIAMETER = "MD"
    _CIRCUMFERENCE = "TC"
    _CONVEX_HULL = "CH"
    _DENSITY = "DE"
    _LABELS = {_AREA: "Total surface area",
               _DIAMETER: "Mean dot diameter",
               _CIRCUMFERENCE: "Total circumference",
               _CONVEX_HULL: "Convex hull",
               _DENSITY: "Density"}

    def __init__(self, type_code, value):
        self._type = type_code
        self.value = value

    @property
    def type_code(self):
        return self._type

    @property
    def long_label(self):
        if self._type is None:
            return "undefined"
        else:
            return _ContinuousProperty._LABELS[self._type]

    def __str__(self):
        return self.long_label

    def __lt__(self, other):
        return self.value < other.value

    def __eq__(self, other):
        return self.value == other.value

    @property
    def _dependencies(self): #todo should be defined in childen classes
        if isinstance(self, (MeanDotDiameter, TotalCircumference, SurfaceArea)):
            return [_ContinuousProperty._DIAMETER,
                    _ContinuousProperty._CIRCUMFERENCE,
                    _ContinuousProperty._AREA]
        elif isinstance(self, ConvexHull):
            return [_ContinuousProperty._CONVEX_HULL]
        elif isinstance(self, Density):
            dep = [_ContinuousProperty._DENSITY]
            if self.match_ratio_convhull2area <1: #area involved
                dep.extend(SurfaceArea()._dependencies)
            if self.match_ratio_convhull2area>0: # convex hull involved
                dep.extend(ConvexHull()._dependencies)
            return dep
        return []


    def is_dependent(self, other):

        return other._type in self._dependencies


class SurfaceArea(_ContinuousProperty):
    """"""
    def __init__(self, value=None):
        _ContinuousProperty.__init__(self, _ContinuousProperty._AREA, value)

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_total_surface_area


class MeanDotDiameter(_ContinuousProperty):
    """"""
    def __init__(self, value=None):
        _ContinuousProperty.__init__(self, _ContinuousProperty._DIAMETER, value)

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_mean_dot_diameter


class TotalCircumference(_ContinuousProperty):
    """"""
    def __init__(self, value=None):
        _ContinuousProperty.__init__(self, _ContinuousProperty._CIRCUMFERENCE, value)

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_total_circumference


class ConvexHull(_ContinuousProperty):
    """"""
    def __init__(self, value=None, match_presision=0.01):
        _ContinuousProperty.__init__(self, _ContinuousProperty._CONVEX_HULL, value)
        self.match_presision = match_presision

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_convex_hull_area


class Density(_ContinuousProperty):
    """"""
    def __init__(self, value=None, match_ratio_convhull2area=0.5, convex_hull_precision=None):

        _ContinuousProperty.__init__(self, _ContinuousProperty._DENSITY, value)
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
        if v<0 or v>1:
            raise RuntimeError("Match_ratio_convhull2area has to be between 0 and 1.")

    def set_value(self, reference_dot_array):
        self.value = reference_dot_array.prop_density


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
            raise TypeError("Parameter is not a continuous properties or a " +\
                            "list of continuous properties")
        elif check_set_value and x.value is None:
            raise RuntimeError("Value of continuous property {} is not defined.".format(
                        type(x).__name__))

    # check dependencies
    for a,b in itertools.combinations(lcp, 2):
        if not a.is_independent(b):
            raise RuntimeError("Incompatible properties to match: {} & {}".format(
                                    type(a).__name__, type(b).__name__))

