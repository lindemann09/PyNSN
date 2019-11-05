"""
Dot Array
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from collections import OrderedDict
from hashlib import md5
import json

import numpy as np
from scipy import spatial
from .misc import log2, numpy_vector, EfficientConvexHullDots
from . import visual_features as vf


TWO_PI = 2 * np.pi

class DotCollection(object):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    def __init__(self, xy=None, diameters=None):

        self._xy = np.array([])
        self._diameters = np.array([])
        self._ch = None
        if (xy, diameters) != (None, None):
            self.append(xy, diameters)
        self.set_array_modified()

    @property
    def xy(self):
        return self._xy

    @property
    def rounded_xy(self):
        """rounded to integer"""
        return np.array(np.round(self._xy), dtype=np.int32)

    @property
    def center_of_outer_positions(self):
        minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def feature_mean_item_diameter(self):
        return np.mean(self._diameters)

    @property
    def feature_total_surface_area(self):
        return np.sum(self.surface_areas)

    @property
    def convex_hull_positions(self):  # FIXME not defined for l<3
        return self._xy[self.convex_hull_indices, :]

    @property
    def convex_hull_indices(self):  # FIXME not defined for l<3
        """this convex_hull takes into account the dot diameter"""
        return self._ch.convex_hull.vertices

    @property
    def surface_areas(self):
        return np.pi * (self._diameters ** 2) / 4.0

    @property
    def perimeter(self):
        return np.pi * self._diameters

    @property
    def feature_mean_item_surface_area(self):
        return np.mean(self.surface_areas)

    @property
    def feature_total_perimeter(self):
        return np.sum(self.perimeter)

    @property
    def feature_mean_item_perimeter(self):
        return np.mean(self.perimeter)

    @property
    def feature_field_area(self):  # todo: not defined for small n
        return self._ch.convex_hull.volume

    @property
    def feature_numerosity(self):
        return len(self._xy)

    @property
    def feature_converage(self):
        """ percent coverage in the field area. It takes thus the item size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity

        """

        try:
            return self.feature_total_surface_area / self.feature_field_area
        except:
            return None

    @property
    def feature_logSize(self):
        return log2(self.feature_total_surface_area) + log2(self.feature_mean_item_surface_area)

    @property
    def feature_logSpacing(self):
        return log2(self.feature_field_area) + log2(self.feature_sparsity)

    @property
    def feature_sparsity(self):
        return self.feature_field_area / self.feature_numerosity

    @property
    def convex_hull_positions_full(self):  # FIXME not defined for l<3
        """this convex_hull takes into account the dot diameter"""
        return self._ch.full_xy

    @property
    def feature_field_area_full(self):  # FIXME not used (correct?)
        return self._ch.full_field_area

    @property
    def object_id(self):
        """md5_hash of position, diameter"""

        m = md5()
        m.update(self._xy.tobytes())  # to byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        m.update(self.surface_areas.tobytes())
        return m.hexdigest()

    def as_dict(self, rounded_values=False):
        if rounded_values:
            xy = np.round(self._xy).astype(np.int).tolist()
            dia = np.round(self._diameters).astype(np.int).tolist()
        else:
            xy = self._xy.tolist()
            dia = self._diameters.tolist()
        return {"object_id": self.object_id, "xy": xy, "diameters": dia}

    def save(self, json_filename, rounded_values=False, indent=None):

        with open(json_filename, 'w') as fl:
            json.dump(self.as_dict(rounded_values), fl, indent=indent)

    def load(self, json_filename):

        with open(json_filename, 'r') as fl:
            dict = json.load(fl)
        self.read_from_dict(dict)

    def read_from_dict(self, dict):
        """read Dot collection from dict"""
        self._xy = np.array(dict["xy"])
        self._diameters = np.array(dict["diameters"])
        self.set_array_modified()

    def get_features_dict(self):
        """ordered dictionary with the most important feature"""
        rtn = [("Object_id", self.object_id),
               ("Numerosity", self.feature_numerosity),
               (vf.TotalSurfaceArea.label, self.feature_total_surface_area),
               (vf.ItemSurfaceArea.label, self.feature_mean_item_surface_area),
               (vf.ItemDiameter.label, self.feature_mean_item_diameter),
               (vf.ItemPerimeter.label, self.feature_mean_item_perimeter),
               (vf.TotalPerimeter.label, self.feature_total_perimeter),
               (vf.FieldArea.label, self.feature_field_area),
               (vf.Sparsity.label, self.feature_sparsity),
               (vf.Coverage.label, self.feature_converage),
               (vf.LogSize.label, self.feature_logSize),
               (vf.LogSpacing.label, self.feature_logSpacing)]
        return OrderedDict(rtn)

    def get_features_text(self, with_object_id=True, extended_format=False, spacing_char="."):
        if extended_format:
            rtn = None
            for k, v in self.get_features_dict().items():
                if rtn is None:
                    if with_object_id:
                        rtn = "- {}\n".format(v)
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
            if with_object_id:
                rtn = "id: {}".format(self.object_id)
            else:
                rtn = ""
            rtn += "n: {}, TSA: {}, ISA: {}, FA: {}, SPAR: {:.3f}, logSIZE: {:.2f}, logSPACE: {:.2f} COV: {:.2f}".format(
                self.feature_numerosity,
                int(self.feature_total_surface_area),
                int(self.feature_mean_item_surface_area),
                int(self.feature_field_area),
                self.feature_sparsity,
                self.feature_logSize,
                self.feature_logSpacing,
                self.feature_converage)
        return rtn

    def set_array_modified(self):
        self._ch = EfficientConvexHullDots(self._xy, self._diameters)

    @property
    def diameters(self):
        return self._diameters

    def append(self, xy, item_diameters):
        """append dots using numpy array"""

        # ensure numpy array
        item_diameters = numpy_vector(item_diameters)
        # ensure xy is a 2d array
        xy = np.array(xy)
        if xy.ndim == 1 and len(xy) == 2:
            xy = xy.reshape((1, 2))
        if xy.ndim != 2:
            raise RuntimeError("Bad shaped data: xy must be pair of xy-values or a list of xy-values")

        if xy.shape[0] != len(item_diameters):
            raise RuntimeError("Bad shaped data: " + u"xy has not the same length as item_diameters")

        if len(self._xy) == 0:
            self._xy = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
        self._xy = np.append(self._xy, xy, axis=0)
        self._diameters = np.append(self._diameters, item_diameters)
        self.set_array_modified()

    def clear(self):
        self._xy = np.array([[]])
        self._diameters = np.array([])
        self.set_array_modified()

    def delete(self, index):
        self._xy = np.delete(self._xy, index, axis=0)
        self._diameters = np.delete(self._diameters, index)
        self.set_array_modified()

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """
        if indices is None:
            indices = list(range(self.feature_numerosity))

        return DotCollection(xy=self._xy[indices, :].copy(),
                             diameters=self._diameters[indices].copy())

    def distances(self, xy, diameter):
        """Distances toward a single point (xy, diameter) """
        if len(self._xy) == 0:
            return np.array([])
        else:
            return np.hypot(self._xy[:, 0] - xy[0], self._xy[:, 1] - xy[1]) - \
                   ((self._diameters + diameter) / 2.0)

    @property
    def center_of_mass(self):
        weighted_sum = np.sum(self._xy * self._diameters[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self._diameters)

    @property
    def feature_mean_item_diameter(self):
        return np.mean(self._diameters)

    def get_distance_matrix(self, between_positions=False):
        """between position ignores the dot size"""
        dist = spatial.distance.cdist(self._xy, self._xy)  # matrix with all distance between all points
        if not between_positions:
            # subtract dot diameter
            radii_mtx = np.ones((self.feature_numerosity, 1)) + self._diameters[:, np.newaxis].T / 2
            dist -= radii_mtx  # for each row
            dist -= radii_mtx.T  # each each column
        return dist

    @property
    def expension(self):
        """ maximal distance between to points plus diameter of the two points"""

        dist = self.get_distance_matrix(between_positions=True)
        # add dot diameter
        radii_mtx = np.ones((self.feature_numerosity, 1)) + self._diameters[:, np.newaxis].T / 2
        dist += radii_mtx  # add to each row
        dist += radii_mtx.T  # add two each column
        return np.max(dist)

    def center_array(self):
        self._xy -= self.center_of_outer_positions
        self.set_array_modified()
