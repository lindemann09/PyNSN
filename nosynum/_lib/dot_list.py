"""
Dot Array
"""
from __future__ import print_function, division, unicode_literals
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import random
from collections import namedtuple
from hashlib import md5
import numpy as np
from scipy.spatial import ConvexHull, distance

TWO_PI = 2 * np.pi

class DotListProperties(namedtuple("DAProperties",["object_id", "numerosity", "mean_dot_diameter",
                                                "total_surface_area", "convex_hull_area",
                                                "density", "total_circumference"])):
    __slots__ = ()

    @classmethod
    def _make_arrays(cls):
        """# create DotListProperties with empty arrays """
        n_prop = len(cls._fields)
        return cls._make(map(lambda _: [], range(n_prop)))

    @classmethod
    def get_np_array_column_names(cls):
        return cls._fields[1:]

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


class DotList(object):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    OBJECT_ID_LENGTH = 8

    def __init__(self, xy_positions=(), diameters=()):

        self.xy = np.array(xy_positions)
        self.diameters = np.array(diameters)

    def clear(self):
        self.xy = np.array([[]])
        self.diameters = np.array([])

    @property
    def object_id(self):
        """md5_hash of position, diameter"""

        m = md5()
        m.update(self.xy.tobytes())
        m.update(self.diameters.tobytes())
        return m.hexdigest()[:DotList.OBJECT_ID_LENGTH]


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
        if len(self.xy) == 0:
            return np.array([])
        else:
            return np.hypot(self.xy[:, 0] - xy[0], self.xy[:, 1] - xy[1]) - \
                   ((self.diameters + diameter) / 2.0)

    @property
    def rounded_xy(self):
        """rounded to integer"""
        return np.array(np.round(self.xy), dtype=np.int32)

    @property
    def rounded_diameters(self):
        """rounded to integer"""
        return np.array(np.round(self.diameters), dtype=np.int32)

    @property
    def center_of_outer_positions(self):
        minmax = np.array((np.min(self.xy, axis=0), np.max(self.xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def center_of_mass(self):
        weighted_sum = np.sum(self.xy * self.diameters[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self.diameters)

    @property
    def surface_areas(self):
        return np.pi * (self.diameters ** 2) / 4.0

    @property
    def circumferences(self):
        return np.pi * self.diameters

    @property
    def convex_hull_positions(self):
        x = ConvexHull(self.xy).vertices
        return self.xy[x, :]

    @property
    def convex_hull(self):
        """this convex_hull takes into account the dot diameter"""
        idx = ConvexHull(self.xy).vertices
        center = self.center_of_outer_positions
        polar_centered = DotList._cartesian2polar(self.xy[idx, :] - center)
        polar_centered[:, 0] = polar_centered[:, 0] + (self.diameters[idx] / 2)
        xy = DotList._polar2cartesian(polar_centered) + center
        return xy

    def _jitter_identical_positions(self, jitter_size=0.1):
        """jitters points with identical position"""

        for idx, ref_dot in enumerate(self.xy):
            identical = np.where(np.all(np.equal(self.xy, ref_dot), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != idx:
                        self.xy[x, :] -= DotList._polar2cartesian([[jitter_size, random.random() * TWO_PI]])[0]

    def _remove_overlap_for_dot(self, dot_id, minimum_gap):
        """remove overlap for one point
        helper function, please use realign
        """

        dist = self.distances(self.xy[dot_id, :], self.diameters[dot_id])

        shift_required = False
        idx = np.where(dist < minimum_gap)[0].tolist()  # overlapping dot ids
        if len(idx) > 1:
            idx.remove(dot_id)  # don't move yourself
            if np.sum(
                    np.all(self.xy[idx,] == self.xy[dot_id, :], axis=1)) > 0:  # check if there is an identical position
                self._jitter_identical_positions()

            tmp_polar = DotList._cartesian2polar(self.xy[idx, :] - self.xy[dot_id, :])
            tmp_polar[:, 0] = 0.000000001 + minimum_gap - dist[idx]  # determine movement size
            xy = DotList._polar2cartesian(tmp_polar)
            self.xy[idx, :] = np.array([self.xy[idx, 0] + xy[:, 0], self.xy[idx, 1] + xy[:, 1]]).T
            shift_required = True

        return shift_required

    @property
    def prop_mean_dot_diameter(self):
        return self.diameters.mean()

    @property
    def prop_total_surface_area(self):
        return np.sum(self.surface_areas)

    @property
    def prop_total_circumference(self):
        return np.sum(self.circumferences)

    @property
    def prop_convex_hull_area_positions(self):
        return ConvexHull(self.xy).area

    @property
    def prop_convex_hull_area(self):
        return ConvexHull(self.convex_hull).area

    @property
    def prop_numerosity(self):
        return len(self.xy)

    @property
    def prop_density(self):
        """density takes into account the convex hull"""
        try:
            return self.prop_convex_hull_area / self.prop_total_surface_area  # todo: positions conved hull
        except:
            return None

    def get_properties(self):
        """returns a named tuple """
        return DotListProperties(self.object_id, self.prop_numerosity, self.prop_mean_dot_diameter,
                                 self.prop_total_surface_area, self.prop_convex_hull_area, self.prop_density,
                                 self.prop_total_circumference)


    def get_distance_matrix(self, between_positions=False):
        """between position ignores the dot size"""
        dist = distance.cdist(self.xy, self.xy)  # matrix with all distance between all points
        if not between_positions:
            # subtract dot diameter
            radii_mtx = np.ones((self.prop_numerosity, 1)) + self.diameters[:, np.newaxis].T / 2
            dist -= radii_mtx  # for each row
            dist -= radii_mtx.T  # each each column
        return dist

    @property
    def expension(self):
        """ maximal distance between to points plus diameter of the two points"""

        dist = self.get_distance_matrix(between_positions=True)
        # add dot diameter
        radii_mtx = np.ones((self.prop_numerosity, 1)) + self.diameters[:, np.newaxis].T / 2
        dist += radii_mtx  # add to each row
        dist += radii_mtx.T  # add two each column
        return np.max(dist)