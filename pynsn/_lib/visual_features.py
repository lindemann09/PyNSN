# calculates visual features of a dot array/ dot cloud

from collections import OrderedDict

import numpy as np
from scipy import spatial
from . import misc
from .geometry import cartesian2polar, polar2cartesian
from . import arrays


class VisualFeatures(object):

    LOG_SIZE = "Log Size"
    TOTAL_SURFACE_AREA = "Total surface area"
    ITEM_DIAMETER = "Mean item diameter"
    ITEM_SURFACE_AREA = "Mean item surface area"
    ITEM_PERIMETER = "Total perimeter"
    TOTAL_PERIMETER = "Mean item perimeter"
    RECT_SIZE = "Mean Rectangle Size"
    LOG_SPACING = "Log Spacing"
    SPARSITY = "Sparsity"
    FIELD_AREA = "Field area"
    COVERAGE = "Coverage"

    SIZE_FEATURES = (LOG_SIZE, TOTAL_SURFACE_AREA, ITEM_DIAMETER,
                     ITEM_SURFACE_AREA, ITEM_PERIMETER, TOTAL_PERIMETER)

    SPACE_FEATURES = (LOG_SPACING, SPARSITY, FIELD_AREA)

    ALL_FEATURES = SIZE_FEATURES + SPACE_FEATURES + (COVERAGE,)

    @staticmethod
    def are_dependent(featureA, featureB):
        """returns true if both features are not independent"""
        for l in [VisualFeatures.SIZE_FEATURES, VisualFeatures.SPACE_FEATURES]:
            if featureA in l and featureB in l:
                return True
        return False

    def __init__(self, object_array):
        # _lib or dot_cloud
        assert isinstance(object_array, (arrays.RectangleArray, arrays.DotArray))
        self.oa = object_array
        self._convex_hull = None

    def reset(self):
        """reset to enforce recalculation of certain parameter """
        self._convex_hull = None

    @property
    def convex_hull(self):
        if self._convex_hull is None:
            self._convex_hull = ConvexHull(self.oa)
        return self._convex_hull

    @property
    def mean_item_diameter(self):
        if not isinstance(self.oa, arrays.DotArray):
            return None
        return np.mean(self.oa.diameters)

    @property
    def mean_rectangle_size(self):
        if not isinstance(self.oa, arrays.RectangleArray):
            return None
        return np.mean(self.oa.sizes, axis=0)


    @property
    def total_surface_area(self):
        return np.sum(self.oa.surface_areas)

    @property
    def mean_item_surface_area(self):
        return np.mean(self.oa.surface_areas)

    @property
    def total_perimeter(self):
        return np.sum(self.oa.perimeter)

    @property
    def mean_item_perimeter(self):
        return np.mean(self.oa.perimeter)

    @property
    def field_area(self):
        return self.convex_hull.convex_hull_position.volume

    @property
    def numerosity(self):
        return len(self.oa._xy)

    @property
    def converage(self):
        """ percent coverage in the field area. It takes thus the item size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity

        """

        try:
            return self.total_surface_area / self.field_area
        except:
            return None

    @property
    def logSize(self):
        return misc.log2(self.total_surface_area) + misc.log2(
            self.mean_item_surface_area)

    @property
    def logSpacing(self):
        return misc.log2(self.field_area) + misc.log2(self.sparsity)

    @property
    def sparsity(self):
        return self.field_area / self.numerosity

    @property
    def field_area_full(self):  # TODO not tested
        return self.convex_hull.outer_field_area

    def get(self, feature):
        """returns a feature"""

       # Adapt
        if feature == VisualFeatures.ITEM_DIAMETER:
            return self.mean_item_diameter

        elif feature == VisualFeatures.ITEM_PERIMETER:
            return self.mean_item_perimeter

        elif feature == VisualFeatures.TOTAL_PERIMETER:
            return self.total_perimeter

        elif feature == VisualFeatures.ITEM_SURFACE_AREA:
            return self.mean_item_surface_area

        elif feature == VisualFeatures.TOTAL_SURFACE_AREA:
            return self.total_surface_area

        elif feature == VisualFeatures.LOG_SIZE:
            return self.logSize

        elif feature == VisualFeatures.LOG_SPACING:
            return self.logSpacing

        elif feature == VisualFeatures.SPARSITY:
            return self.sparsity

        elif feature == VisualFeatures.FIELD_AREA:
            return self.field_area

        elif feature == VisualFeatures.COVERAGE:
            return self.converage

        else:
            raise ValueError("{} is a unkown visual feature".format(feature))


    def as_dict(self):
        """ordered dictionary with the most important feature"""
        rtn = [("Hash", self.oa.hash),
               ("Numerosity", self.numerosity),
               (VisualFeatures.TOTAL_SURFACE_AREA, self.total_surface_area),
               (VisualFeatures.ITEM_SURFACE_AREA, self.mean_item_surface_area),
               ("?", None), # placeholder
               (VisualFeatures.ITEM_PERIMETER, self.mean_item_perimeter),
               (VisualFeatures.TOTAL_PERIMETER, self.total_perimeter),
               (VisualFeatures.FIELD_AREA, self.field_area),
               (VisualFeatures.SPARSITY, self.sparsity),
               (VisualFeatures.COVERAGE, self.converage),
               (VisualFeatures.LOG_SIZE, self.logSize),
               (VisualFeatures.LOG_SPACING, self.logSpacing)]

        if isinstance(self.oa, arrays.DotArray):
            rtn[4] = (VisualFeatures.ITEM_DIAMETER, self.mean_item_diameter)
        elif isinstance(self.oa, arrays.RectangleArray):
            rtn[4] = (VisualFeatures.RECT_SIZE, self.mean_rectangle_size)
        return OrderedDict(rtn)

    def __str__(self):
        return self.as_text()

    def as_text(self, with_hash=True, extended_format=False, spacing_char="."):
        if extended_format:
            rtn = None
            for k, v in self.as_dict().items():
                if rtn is None:
                    if with_hash:
                        rtn = "- {}: {}\n".format(k, v)
                    else:
                        rtn = ""
                else:
                    if rtn == "":
                        name = "- " + k
                    else:
                        name = "  " + k
                    try:
                        value = "{0:.2f}\n".format(v)  # try rounding
                    except:
                        value = "{}\n".format(v)

                    rtn += name + (spacing_char * (22 - len(name))) + (" " * (14 - len(value))) + value
        else:
            if with_hash:
                rtn = "ID: {} ".format(self.oa.hash)
            else:
                rtn = ""
            rtn += "N: {}, TSA: {}, ISA: {}, FA: {}, SPAR: {:.3f}, logSIZE: " \
                   "{:.2f}, logSPACE: {:.2f} COV: {:.2f}".format(
                self.numerosity,
                int(self.total_surface_area),
                int(self.mean_item_surface_area),
                int(self.field_area),
                self.sparsity,
                self.logSize,
                self.logSpacing,
                self.converage)
        return rtn


class ConvexHull(object):
    """convenient wrapper class for calculation of convex hulls
    """

    def __init__(self, object_array):
        assert isinstance(object_array, (arrays.RectangleArray, arrays.DotArray))

        self._xy = object_array.xy

        if isinstance(object_array, arrays.DotArray):
            # centered polar coordinates
            minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)  # center outer positions
            xy_centered = self._xy - center
            # outer positions
            polar_centered = cartesian2polar(xy_centered)
            polar_centered[:, 0] = polar_centered[:, 0] + (object_array.diameters / 2)
            self._outer_xy = polar2cartesian(polar_centered) + center

        elif isinstance(object_array, arrays.RectangleArray):
            # get all edges
            edges = []
            for rect in object_array.get():
                edges.extend([e.xy for e in rect.edges()])
            self._outer_xy = np.array(edges)

        self.convex_hull_position = spatial.ConvexHull(self._xy)
        self.convex_hull_outer = spatial.ConvexHull(self._outer_xy)

    @property
    def position_xy(self):
        return self._xy[self.convex_hull_position.vertices, :]

    @property
    def outer_xy(self):
        return self._outer_xy[self.convex_hull_outer.vertices, :]

    @property
    def outer_field_area(self):
        return self.convex_hull_outer.volume

    @property
    def position_field_area(self):
        return self.convex_hull_position.volume

