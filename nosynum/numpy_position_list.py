"""
Dot Array
"""
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'


import numpy as np
from .dot import Dot

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
