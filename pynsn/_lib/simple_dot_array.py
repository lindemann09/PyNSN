"""
Dot Array
"""
from __future__ import print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import random
from copy import copy
from collections import namedtuple
from hashlib import md5
import numpy as np
from scipy import spatial
from .item_features import numpy_vector
from . import continuous_property as cp

TWO_PI = 2 * np.pi

class DotArrayProperties(namedtuple("DAProperties", ["object_id", "numerosity", "mean_dot_diameter",
                                                "total_surface_area", "convex_hull_area",
                                                "density", "total_circumference"])):
    __slots__ = ()

    @classmethod
    def _make_arrays(cls):
        """# create DotListProperties with empty arrays """
        n_prop = len(cls._fields)
        return cls._make(map(lambda _: [], range(n_prop)))

    @property
    def property_names(self):
        return self._fields[1:]

    @property
    def np_array(self):
        return np.array(self[1:]).T

    @property
    def size(self):
        try:
            return len(self.numerosity)
        except:
            return 1

    def get_csv(self, variable_names=True):
        rtn = ""
        if variable_names:
            rtn += u", ".join(self._fields) + "\n"
        if self.size>1:
            for i, prop in enumerate(self.np_array):
                rtn += self.object_id[i] + u", " + ", ".join(map(str, prop)) + "\n"
        else:
            rtn += u", ".join(map(str, self)) + u"\n"

        return rtn

    @property
    def split(self):
        """list of DotArrayProperties objects instead
        one DotArrayProperties object with lists for each property"""
        rtn = []
        for i in range(len(self[0])):
            tmp = []
            for f in range(len(self._fields)):
                tmp.append(self[f][i])
            rtn.append(DotArrayProperties._make(tmp))
        return rtn

    def get_nice_text(self, spacing_char="."):
        txt = None
        for k,v in self._asdict().items():
            if txt is None:
                txt = "- {}\n".format(self.object_id)
            else:
                name = "  " + k
                try:
                    value = "{0:.2f}\n".format(v) # try rounding
                except:
                    value = "{}\n".format(v)

                txt += name + (spacing_char*(23-len(name))) + (" "*(10-len(value))) + value
        return txt

class SimpleDotArray(object):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    OBJECT_ID_LENGTH = 8

    def __init__(self, xy=None, diameters=None):

        self._xy = np.array([])
        self._diameters = np.array([])
        self._ch = EfficientConvexHull(self._xy, self._diameters)
        if (xy, diameters) != (None, None):
            self.append(xy, diameters)

    @property
    def xy(self):
        return self._xy

    @property
    def diameters(self):
        return self._diameters

    def append(self, xy, diameters):
        """append dots using numpy array"""

        # ensure numpy array
        diameters = numpy_vector(diameters)
        # ensure xy is a 2d array
        xy = np.array(xy)
        if xy.ndim == 1 and len(xy) == 2:
            xy = xy.reshape((1, 2))
        if xy.ndim != 2:
            raise RuntimeError("Bad shaped data: xy must be pair of xy-values or a list of xy-values")

        if xy.shape[0] != len(diameters):
            raise RuntimeError("Bad shaped data: " + u"xy has not the same length as diameter")

        if len(self._xy) == 0:
            self._xy = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
        self._xy = np.append(self._xy, xy, axis=0)
        self._diameters = np.append(self._diameters, diameters)
        self._ch = EfficientConvexHull(self._xy, self._diameters)

    def clear(self):
        self._xy = np.array([[]])
        self._diameters = np.array([])
        self._ch = EfficientConvexHull(self._xy, self._diameters)

    def delete(self, index):
        self._xy = np.delete(self._xy, index, axis=0)
        self._diameters = np.delete(self._diameters, index)
        self._ch = EfficientConvexHull(self._xy, self._diameters)

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """
        if indices is None:
            indices = list(range(self.prop_numerosity))

        return SimpleDotArray(xy=self._xy[indices, :].copy(),
                              diameters=self._diameters[indices].copy())


    @property
    def object_id(self):
        """md5_hash of position, diameter"""

        m = md5()
        m.update(self._xy.tobytes())
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
    def circumferences(self):
        return np.pi * self._diameters

    @property
    def convex_hull_positions(self):
        x = self._ch.convex_hull_object.vertices
        return self._xy[x, :]

    @property
    def convex_hull(self): #FIXME
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
                    np.all(self._xy[idx,] == self._xy[dot_id, :], axis=1)) > 0:  # check if there is an identical position
                self._jitter_identical_positions()

            tmp_polar = SimpleDotArray._cartesian2polar(self._xy[idx, :] - self._xy[dot_id, :])
            tmp_polar[:, 0] = 0.000000001 + minimum_gap - dist[idx]  # determine movement size
            xy = SimpleDotArray._polar2cartesian(tmp_polar)
            self._xy[idx, :] = np.array([self._xy[idx, 0] + xy[:, 0], self._xy[idx, 1] + xy[:, 1]]).T
            shift_required = True

        return shift_required

    @property
    def prop_mean_dot_diameter(self):
        return self._diameters.mean()

    @property
    def prop_total_surface_area(self):
        return np.sum(self.surface_areas)

    @property
    def prop_total_circumference(self):
        return np.sum(self.circumferences)

    @property
    def prop_convex_hull_area_positions(self):
        return self._ch.convex_hull_object.area

    @property
    def prop_convex_hull_area(self):
        return self._ch.full_area

    @property
    def prop_numerosity(self):
        return len(self._xy)

    @property
    def prop_density(self):
        """density takes into account the convex hull"""
        try:
            return self.prop_convex_hull_area / self.prop_total_surface_area  # todo: positions conved hull
        except:
            return None

    def get_properties(self):
        """returns a named tuple """
        return DotArrayProperties(self.object_id, self.prop_numerosity, self.prop_mean_dot_diameter,
                                  self.prop_total_surface_area, self.prop_convex_hull_area, self.prop_density,
                                  self.prop_total_circumference)


    def get_distance_matrix(self, between_positions=False):
        """between position ignores the dot size"""
        dist = spatial.distance.cdist(self._xy, self._xy)  # matrix with all distance between all points
        if not between_positions:
            # subtract dot diameter
            radii_mtx = np.ones((self.prop_numerosity, 1)) + self._diameters[:, np.newaxis].T / 2
            dist -= radii_mtx  # for each row
            dist -= radii_mtx.T  # each each column
        return dist

    @property
    def expension(self):
        """ maximal distance between to points plus diameter of the two points"""

        dist = self.get_distance_matrix(between_positions=True)
        # add dot diameter
        radii_mtx = np.ones((self.prop_numerosity, 1)) + self._diameters[:, np.newaxis].T / 2
        dist += radii_mtx  # add to each row
        dist += radii_mtx.T  # add two each column
        return np.max(dist)


    def match(self, match_properties, center_array=True, match_dot_array=None):
        """
        match_properties: continuous property or list of continuous properties
        several properties to be matched

        if match dot array is specified, array will be match to match_dot_array, otherwise
        the values defined in match_properties are used.
        """

        # type check
        if not isinstance(match_properties, (list, tuple)):
            match_properties = [match_properties]

        cp.check_list_continuous_properties(match_properties, check_set_value=match_dot_array is None)

        # copy and change values to match this stimulus
        if match_dot_array is None:
            match_props = match_properties
        else:
            match_props = []
            for m in match_properties:
                m = copy(m)
                m.set_value(match_dot_array)
                match_props.append(m)

        # Adapt
        for cont_prop in match_props:
            if isinstance(cont_prop, cp.SurfaceArea):
                self._match_total_surface_area(surface_area=cont_prop.value)
            elif isinstance(cont_prop, cp.DotDiameter):
                self._match_mean_dot_diameter(mean_dot_diameter=cont_prop.value)
            elif isinstance(cont_prop, cp.Circumference):
                self._match_total_circumference(total_circumference=cont_prop.value)
            elif isinstance(cont_prop, cp.ConvexHull):
                self._match_convex_hull_area(convex_hull_area=cont_prop.value,
                                             precision=cont_prop.match_presision)
            elif isinstance(cont_prop, cp.Density):
                self._match_density(density=cont_prop.value,
                                    precision=cont_prop.convex_hull_precision,
                                    adaptation_CH2TA_ratio=cont_prop.match_ratio_convhull2area)

        if center_array:
            self._xy -= self.center_of_outer_positions

    def _match_total_surface_area(self, surface_area):
        # changes diameter
        a_scale = (surface_area / self.prop_total_surface_area)
        self._diameters = np.sqrt(self.surface_areas * a_scale) * 2 / np.sqrt(
            np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)

    def _match_mean_dot_diameter(self, mean_dot_diameter):
        # changes diameter

        scale = mean_dot_diameter / self.prop_mean_dot_diameter
        self._diameters = self._diameters * scale

    def _match_total_circumference(self, total_circumference):
        # linear to fit_mean_dot_diameter, but depends on numerosity
        mean_dot_diameter = total_circumference / (self.prop_numerosity * np.pi)
        self._match_mean_dot_diameter(mean_dot_diameter)

    def _match_convex_hull_area(self, convex_hull_area, precision=0.001):
        """changes the convex hull area to a desired size with certain precision

        iterative method can takes some time.
        """

        current = self.prop_convex_hull_area

        if current is None:
            return  # not defined

        # iteratively determine scale
        scale = 1  # find good scale
        step = 0.1
        if convex_hull_area < current:  # current too larger
            step *= -1

        # centered  points
        old_center = self.center_of_outer_positions
        centered_polar = SimpleDotArray._cartesian2polar(self._xy - old_center)

        while abs(current - convex_hull_area) > precision:
            scale += step

            self._xy = SimpleDotArray._polar2cartesian(centered_polar * [scale, 1])
            current = self.prop_convex_hull_area

            if (current < convex_hull_area and step < 0) or \
                    (current > convex_hull_area and step > 0):
                step *= -0.2  # change direction and finer grain

        self._xy += old_center

    def _match_density(self, density, precision=0.001, adaptation_CH2TA_ratio=0.5):
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """

        # dens = convex_hull_area / total_surface_area
        if adaptation_CH2TA_ratio < 0 or adaptation_CH2TA_ratio > 1:
            adaptation_CH2TA_ratio = 0.5

        area_change100 = (self.prop_convex_hull_area / density) - self.prop_total_surface_area
        d_change_area = area_change100 * (1 - adaptation_CH2TA_ratio)
        if abs(d_change_area) > 0:
            self._match_total_surface_area(surface_area=self.prop_total_surface_area + d_change_area)

        self._match_convex_hull_area(convex_hull_area=self.prop_total_surface_area * density,
                                     precision=precision)


    def remove_overlap_from_inner_to_outer(self, minimum_gap):

        shift_required = False
        # from inner to outer remove overlaps
        for i in np.argsort(SimpleDotArray._cartesian2polar(self._xy, radii_only=True)):
            if self._remove_overlap_for_dot(dot_id=i, minimum_gap=minimum_gap):
                shift_required = True

        return shift_required

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
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2) # center outer positions

            polar_centered = SimpleDotArray._cartesian2polar(self._xy[idx, :] - center)
            polar_centered[:, 0] = polar_centered[:, 0] + (self._diameters[idx] / 2)
            self._full_xy = SimpleDotArray._polar2cartesian(polar_centered) + center

        return self._full_xy

    @property
    def full_area(self):
        if self._full_area is None:
            self._full_area = spatial.ConvexHull(self.full_xy).area
        return self._full_area