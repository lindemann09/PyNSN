"""
Dot Array
"""
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import random
import numpy as np
from .dot import Dot
TWO_PI = 2*np.pi

class NumpyPositionList(object):
    """Numpy Position list

    for optimized calculations
    """

    def __init__(self, list_of_Dots):
        xy = []
        d = []
        for dot in list_of_Dots:
            xy.append(dot.xy)
            d.append(dot.diameter)

        self._xy = np.array(xy)
        self.diameter = np.array(d)
        self._polar = None

    @property
    def xy(self):
        """numpy array of [x, y]"""

        if self._xy is None:
            self.calc_cartesian()
        return self._xy

    @xy.setter
    def xy(self, xy_position_list):
        self._xy = np.array(xy_position_list)
        self._polar = None

    @property
    def polar(self):
        """numpy array of [radius, angle]"""
        if self._polar is None:
            self.calc_polar()
        return self._polar

    @polar.setter
    def polar(self, radius_angle_list):
        self._polar = np.array(radius_angle_list)
        self._xy = None

    def calc_cartesian(self):
        self._xy = np.array([self._polar[:, 0] * np.cos(self._polar[:, 1]),
                             self._polar[:, 0] * np.sin(self._polar[:, 1])]).T

    def calc_polar(self):
        self._polar = np.array([np.hypot(self._xy[:, 0], self._xy[:, 1]),
                                np.arctan2(self._xy[:, 0], self._xy[:, 1])]).T

    def distance(self, dot):
        return np.hypot(self.xy[:,0] - dot.x, self.xy[:, 1] - dot.y) - \
               ((self.diameter + dot.diameter) / 2.0)


    def jitter_identical_positions(self):
        """jitters points with identical position"""
        tmp = Dot(diameter=0)
        for i in range(len(self._xy)):
            ref_dot = Dot(x=self._xy[i, 0], y=self._xy[i, 1], diameter=self.diameter[i])
            identical = np.where(np.all(np.equal(self._xy, ref_dot.xy), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != i:
                        tmp.polar = (0.1, random.random() * TWO_PI)
                        self._xy[x, :] -= tmp.xy
                        self._polar = None

    def remove_overlap_for_dot(self, dotid, min_gap):
        """remove overlap for one point"""

        ref_dot = Dot(x=self._xy[dotid, 0], y=self._xy[dotid, 1], diameter=self.diameter[dotid])
        dist = self.distance(ref_dot)
        tmp = Dot(diameter=0)
        shift_required = False
        for x in np.where(dist < min_gap)[0]: # x=overlapping dot id
            if x != dotid:
                shift_required = True
                # calc vector (tmp) to shift
                tmp.xy = self._xy[x, :] - ref_dot.xy
                tmp.pos_radius = 0.000000001 + min_gap - dist[x]

                self._xy[x, :] = (self._xy[x, 0] + tmp.x, self._xy[x, 1] + tmp.y)
                shift_required = True

        if shift_required:
            self._polar = None

        return shift_required

