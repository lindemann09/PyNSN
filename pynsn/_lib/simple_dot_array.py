"""
Dot Array
"""
from __future__ import print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import random
from copy import copy
from collections import OrderedDict
from hashlib import md5
import math
import numpy as np
from scipy import spatial
from .item_attributes import numpy_vector
from . import features as cp

TWO_PI = 2 * np.pi

#
# class DotArrayFeature(namedtuple("DAProperties", ["object_id", "numerosity", "item_diameter",
#                                                      "total_surface_area", "field_area",
#                                                      "coverage", "total_perimeter",
#                                                     "sparsity"])): # FIXME use spacity and new feature names
#     __slots__ = ()
#
#     @classmethod
#     def _make_arrays(cls):
#         """# create DotListProperties with empty arrays """
#         n_prop = len(cls._fields)
#         return cls._make(map(lambda _: [], range(n_prop)))
#
#     @property
#     def feature_names(self):
#         return self._fields[1:]
#
#     @property
#     def np_array(self):
#         return np.array(self[1:]).T
#
#     @property
#     def size(self):
#         try:
#             return len(self.numerosity)
#         except:
#             return 1
#
#     def get_csv(self, variable_names=True):
#         rtn = ""
#         if variable_names:
#             rtn += u", ".join(self._fields) + "\n"
#         if self.size > 1:
#             for i, prop in enumerate(self.np_array):
#                 rtn += self.object_id[i] + u", " + ", ".join(map(str, prop)) + "\n"
#         else:
#             rtn += u", ".join(map(str, self)) + u"\n"
#
#         return rtn
#
#     @property
#     def split(self):
#         """list of DotArrayProperties objects instead
#         one DotArrayProperties object with lists for each property"""
#         rtn = []
#         for i in range(len(self[0])):
#             tmp = []
#             for f in range(len(self._fields)):
#                 tmp.append(self[f][i])
#             rtn.append(DotArrayFeature._make(tmp))
#         return rtn

class SimpleDotArray(object):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    OBJECT_ID_LENGTH = 8

    def __init__(self, xy=None, diameters=None):

        self._xy = np.array([])
        self._diameters = np.array([])
        self._ch = None
        if (xy, diameters) != (None, None):
            self.append(xy, diameters)
        self.set_array_modified()

    def set_array_modified(self):
        self._ch = EfficientConvexHull(self._xy, self._diameters)

    @property
    def xy(self):
        return self._xy

    @property
    def item_diameters(self):
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

        return SimpleDotArray(xy=self._xy[indices, :].copy(),
                              diameters=self._diameters[indices].copy())

    @property
    def object_id(self):
        """md5_hash of position, diameter"""

        m = md5()
        m.update(self._xy.tobytes()) # to byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        m.update(self._diameters.tobytes())
        return m.hexdigest()[:SimpleDotArray.OBJECT_ID_LENGTH]

    @staticmethod
    def _polar2cartesian(polar):
        """polar is an 2d-array representing polar coordinates (radius, angle)"""
        polar = np.array(polar)
        return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                         polar[:, 0] * np.sin(polar[:, 1])]).T

    @staticmethod
    def _cartesian2polar(xy, radii_only=False):
        xy = np.array(xy)
        radii = np.hypot(xy[:, 0], xy[:, 1])
        if radii_only:
            return radii
        else:
            return np.array([radii, np.arctan2(xy[:, 1], xy[:, 0])]).T

    def distances(self, xy, diameter):
        """Distances toward a single point (xy, diameter) """
        if len(self._xy) == 0:
            return np.array([])
        else:
            return np.hypot(self._xy[:, 0] - xy[0], self._xy[:, 1] - xy[1]) - \
                   ((self._diameters + diameter) / 2.0)

    @property
    def rounded_xy(self):
        """rounded to integer"""
        return np.array(np.round(self._xy), dtype=np.int32)

    @property
    def rounded_diameters(self):
        """rounded to integer"""
        return np.array(np.round(self._diameters), dtype=np.int32)

    @property
    def center_of_outer_positions(self):
        minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def center_of_mass(self):
        weighted_sum = np.sum(self._xy * self._diameters[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self._diameters)

    @property
    def surface_areas(self):
        return np.pi * (self._diameters ** 2) / 4.0

    @property
    def perimeter(self):
        return np.pi * self._diameters

    @property
    def convex_hull_positions(self):  # FIXME not defined for l<3
        return self._xy[self.convex_hull_indices, :]

    @property
    def convex_hull_indices(self):  # FIXME not defined for l<3
        """this convex_hull takes into account the dot diameter"""
        return self._ch.convex_hull_object.vertices

    @property
    def full_convex_hull_positions(self):  # FIXME not defined for l<3
        """this convex_hull takes into account the dot diameter"""
        return self._ch.full_xy

    def _jitter_identical_positions(self, jitter_size=0.1):
        """jitters points with identical position"""

        for idx, ref_dot in enumerate(self._xy):
            identical = np.where(np.all(np.equal(self._xy, ref_dot), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != idx:
                        self._xy[x, :] -= SimpleDotArray._polar2cartesian([[jitter_size, random.random() * TWO_PI]])[0]

    def _remove_overlap_for_dot(self, dot_id, minimum_gap):
        """remove overlap for one point
        helper function, please use realign
        """

        dist = self.distances(self._xy[dot_id, :], self._diameters[dot_id])

        shift_required = False
        idx = np.where(dist < minimum_gap)[0].tolist()  # overlapping dot ids
        if len(idx) > 1:
            idx.remove(dot_id)  # don't move yourself
            if np.sum(
                    np.all(self._xy[idx,] == self._xy[dot_id, :],
                           axis=1)) > 0:  # check if there is an identical position
                self._jitter_identical_positions()

            tmp_polar = SimpleDotArray._cartesian2polar(self._xy[idx, :] - self._xy[dot_id, :])
            tmp_polar[:, 0] = 0.000000001 + minimum_gap - dist[idx]  # determine movement size
            xy = SimpleDotArray._polar2cartesian(tmp_polar)
            self._xy[idx, :] = np.array([self._xy[idx, 0] + xy[:, 0], self._xy[idx, 1] + xy[:, 1]]).T
            shift_required = True

        return shift_required

    @property
    def feature_item_diameter(self):
        return np.mean(self._diameters)

    @property
    def feature_item_surface_area(self):
        return np.mean(self.surface_areas)

    @property
    def feature_total_surface_area(self):
        return np.sum(self.surface_areas)

    @property
    def feature_total_perimeter(self):
        return np.sum(self.perimeter)

    @property
    def feature_field_area(self): # todo: not defined for small n
        return self._ch.convex_hull_object.volume

    @property
    def feature_field_area_outer(self):  # FIXME not used (correct?)
        return self._ch.full_field_area

    @property
    def feature_numerosity(self):
        return len(self._xy)

    @property
    def feature_converage(self):
        """density takes into account the convex hull"""
        try:
            return self.feature_total_surface_area / self.feature_field_area
        except:
            return None

    @property
    def feature_Size(self):
        return (self.feature_total_surface_area / math.sqrt(self.feature_numerosity)) ** 2

    @property
    def feature_Space(self):
        return (self.feature_field_area / math.sqrt(self.feature_numerosity)) ** 2

    @property
    def feature_sparsity(self):
        return self.feature_field_area / self.feature_numerosity

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

    def match(self, match_features, center_array=True, match_dot_array=None):
        """
        match_properties: continuous property or list of continuous properties
        several properties to be matched

        if match dot array is specified, array will be match to match_dot_array, otherwise
        the values defined in match_properties are used.
        """

        # type check
        if not isinstance(match_features, (list, tuple)):
            match_features = [match_features]

        cp.check_feature_list(match_features, check_set_value=match_dot_array is None)

        # copy and change values to match this stimulus
        if match_dot_array is None:
            match_feat = match_features
        else:
            match_feat = []
            for m in match_features:
                m = copy(m)
                m.set_value(match_dot_array)
                match_feat.append(m)

        # Adapt
        for feat in match_feat:
            if isinstance(feat, cp.ItemDiameter):
                self._match_item_diameter(mean_item_diameter=feat.value)
            elif isinstance(feat, cp.TotalPerimeter):
                mean_dot_diameter = feat.value / (self.feature_numerosity * np.pi)
                self._match_item_diameter(mean_dot_diameter)
            elif isinstance(feat, cp.TotalSurfaceArea):
                self._match_total_surface_area(surface_area=feat.value)
            elif isinstance(feat, cp.FieldArea):
                self._match_field_area(field_area=feat.value,
                                       precision=feat.match_presision)
            elif isinstance(feat, cp.Coverage):
                self._match_coverage(coverage=feat.value,
                                     precision=feat.convex_hull_precision,
                                     match_FA2TA_ratio=feat.match_ratio_fieldarea2totalarea)

        if center_array:
            self._xy -= self.center_of_outer_positions
            self.set_array_modified()

    def _match_total_surface_area(self, surface_area):
        # changes diameter
        a_scale = (surface_area / self.feature_total_surface_area)
        self._diameters = np.sqrt(self.surface_areas * a_scale) * 2 / np.sqrt(
            math.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
        self.set_array_modified()

    def _match_item_diameter(self, mean_item_diameter):
        # changes diameter

        scale = mean_item_diameter / self.feature_item_diameter
        self._diameters = self._diameters * scale
        self.set_array_modified()

    def _match_field_area(self, field_area, precision=0.001):
        """changes the convex hull area to a desired size with certain precision

        iterative method can takes some time.
        """
        current = self.feature_field_area

        if current is None:
            return  # not defined

        # iteratively determine scale
        scale = 1  # find good scale
        step = 0.1
        if field_area < current:  # current too larger
            step *= -1

        # centered points
        old_center = self.center_of_outer_positions
        self._xy -= old_center
        centered_polar = SimpleDotArray._cartesian2polar(self._xy)

        while abs(current - field_area) > precision:

            scale += step

            self._xy = SimpleDotArray._polar2cartesian(centered_polar * [scale, 1])
            self.set_array_modified()  # required to recalc convex hull
            current = self.feature_field_area

            if (current < field_area and step < 0) or \
                    (current > field_area and step > 0):
                step *= -0.2  # change direction and finer grain

        self._xy += old_center
        self.set_array_modified()


    def _match_coverage(self, coverage, precision=0.001, match_FA2TA_ratio=0.5): # FIXME check drifting outwards if extra space is small and match_FA2TA_ratio=1
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """

        # dens = convex_hull_area / total_surface_area
        if match_FA2TA_ratio < 0 or match_FA2TA_ratio > 1:
            match_FA2TA_ratio = 0.5

        total_area_change100 = (coverage * self.feature_field_area ) - self.feature_total_surface_area
        d_change_total_area = total_area_change100 * (1 - match_FA2TA_ratio)
        if abs(d_change_total_area) > 0:
            self._match_total_surface_area(surface_area=self.feature_total_surface_area + d_change_total_area)

        self._match_field_area(field_area=self.feature_total_surface_area / coverage,
                               precision=precision)

    def center_array(self):
        self._xy -= self.center_of_outer_positions
        self.set_array_modified()

    def remove_overlap_from_inner_to_outer(self, minimum_gap):

        shift_required = False
        # from inner to outer remove overlaps
        for i in np.argsort(SimpleDotArray._cartesian2polar(self._xy, radii_only=True)):
            if self._remove_overlap_for_dot(dot_id=i, minimum_gap=minimum_gap):
                shift_required = True

        if shift_required:
            self.set_array_modified()

        return shift_required


    def get_features_dict(self):
        """ordered dictionary with the most important feature"""
        rtn = [("object_id", self.object_id),
               ("numerosity", self.feature_numerosity),
               ("total surface area", self.feature_total_surface_area),
               ("item surface area", self.feature_item_surface_area),
               ("field area", self.feature_field_area),
               ("sparsity", self.feature_sparsity),
               ("converage", self.feature_converage)]
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

                    rtn += name + (spacing_char * (23 - len(name))) + (" " * (10 - len(value))) + value
        else:
            if with_object_id:
                rtn = "id: {}".format(self.object_id)
            else:
                rtn =""
            rtn += "n: {}, TSA: {}, ISA: {}, FA: {}, SPAR: {:.3f}, COV: {:.3f}".format(self.feature_numerosity,
                                                                         int(self.feature_total_surface_area),
                                                                         int(self.feature_item_surface_area),
                                                                         int(self.feature_field_area),
                                                                         self.feature_sparsity,
                                                                         self.feature_converage) # TODO cardinal
        return rtn


########## helper

class EfficientConvexHull(object):
    """helper class for efficent (and not repeated) calulations convex hull"""

    def __init__(self, xy, diameters):
        self._xy = xy
        self._diameters = diameters
        self._positions_ch_object = None
        self._full_xy = None
        self._full_area = None

    @property
    def convex_hull_object(self):
        if self._positions_ch_object is None:
            self._positions_ch_object = spatial.ConvexHull(self._xy)

        return self._positions_ch_object

    @property
    def full_xy(self):
        """this convex_hull takes into account the dot diameter"""
        if self._full_xy is None:
            idx = self.convex_hull_object.vertices

            minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)  # center outer positions

            polar_centered = SimpleDotArray._cartesian2polar(self._xy[idx, :] - center)
            polar_centered[:, 0] = polar_centered[:, 0] + (self._diameters[idx] / 2)
            self._full_xy = SimpleDotArray._polar2cartesian(polar_centered) + center

        return self._full_xy

    @property
    def full_field_area(self):
        if self._full_area is None:
            self._full_area = spatial.ConvexHull(self.full_xy).volume

        return self._full_area
