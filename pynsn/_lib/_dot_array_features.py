# calculates visual features of a dot array/ dot cloud

from collections import OrderedDict

import numpy as np
from scipy import spatial
from . import _misc
from . import features

class DotArrayFeatures(object):

    def __init__(self, dot_array):
        # dot_array or dot_cloud

        self.da = dot_array

    @property
    def mean_item_diameter(self):
        return np.mean(self.da._diameters)

    @property
    def total_surface_area(self):
        return np.sum(self.da.surface_areas)

    @property
    def mean_item_surface_area(self):
        return np.mean(self.da.surface_areas)

    @property
    def total_perimeter(self):
        return np.sum(self.da.perimeter)

    @property
    def mean_item_perimeter(self):
        return np.mean(self.da.perimeter)

    @property
    def field_area(self):
        return self.da.convex_hull.scipy_convex_hull.volume

    @property
    def numerosity(self):
        return len(self.da._xy)

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
        return _misc.log2(self.total_surface_area) + _misc.log2(
            self.mean_item_surface_area)

    @property
    def logSpacing(self):
        return _misc.log2(self.field_area) + _misc.log2(self.sparsity)

    @property
    def sparsity(self):
        return self.field_area / self.numerosity

    @property
    def field_area_full(self):  # FIXME not used (correct?)
        return self.da._ch.full_field_area

    def _get_distance_matrix(self, between_positions=False):
        """between position ignores the dot size"""
        dist = spatial.distance.cdist(self.da._xy, self.da._xy)  #
        # matrix with all distance between all points
        if not between_positions:
            # subtract dot diameter
            radii_mtx = np.ones((self.numerosity, 1)) + \
                        self.da._diameters[:, np.newaxis].T / 2
            dist -= radii_mtx  # for each row
            dist -= radii_mtx.T  # each each column
        return dist

    @property
    def expension(self):
        """ maximal distance between to points plus diameter of the two points"""

        dist = self._get_distance_matrix(between_positions=True)
        # add dot diameter
        radii_mtx = np.ones((self.numerosity, 1)) + self.da._diameters[:,
                                                    np.newaxis].T / 2
        dist += radii_mtx  # add to each row
        dist += radii_mtx.T  # add two each column
        return np.max(dist)

    @property
    def featureXX_coverage_target_area(self):
        """density takes into account the full possible target area (i.e., stimulus radius) """
        try:
            return np.pi * self.da.target_array_radius ** 2 / \
                self.total_surface_area
        except:
            return None # dot defined for DotCloud with no target_array_radius


    def get_features_dict(self):
        """ordered dictionary with the most important feature"""
        rtn = [("Hash", self.da.hash),
               ("Numerosity", self.numerosity),
               (features.TOTAL_SURFACE_AREA, self.total_surface_area),
               (features.ITEM_SURFACE_AREA, self.mean_item_surface_area),
               (features.ITEM_DIAMETER, self.mean_item_diameter),
               (features.ITEM_PERIMETER, self.mean_item_perimeter),
               (features.TOTAL_PERIMETER, self.total_perimeter),
               (features.FIELD_AREA, self.field_area),
               (features.SPARSITY, self.sparsity),
               (features.COVERAGE, self.converage),
               (features.LOG_SIZE, self.logSize),
               (features.LOG_SPACING, self.logSpacing)]
        return OrderedDict(rtn)

    def get_features_text(self, with_hash=True, extended_format=False, spacing_char="."):
        if extended_format:
            rtn = None
            for k, v in self.get_features_dict().items():
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
                rtn = "ID: {} ".format(self.da.hash)
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



