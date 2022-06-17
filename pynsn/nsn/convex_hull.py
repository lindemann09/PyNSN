"""
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from scipy import spatial
from ..lib.geometry import cartesian2polar, polar2cartesian

class ConvexHull(object):
    """convenient wrapper class for calculation of convex hulls
    """

    def __init__(self, xy):
        self._xy = xy
        self.scipy_convex_hull = spatial.ConvexHull(self._xy)

    @property
    def indices(self):
        return self.scipy_convex_hull.vertices

    @property
    def xy(self):
        return self._xy[self.indices, :]


class ConvexHullDots(ConvexHull):
    """convenient wrapper class for calculation of convex hulls
    """

    def __init__(self, xy, diameters):
        ConvexHull.__init__(self, xy=xy)
        self._diameters = diameters
        self._full_xy = None
        self._full_area = None

    @property
    def full_xy(self):
        """this convex_hull takes into account the dot diameter"""
        if self._full_xy is None:
            idx = self.scipy_convex_hull.vertices

            minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)  # center outer positions

            polar_centered = cartesian2polar(self._xy[idx, :] - center)
            polar_centered[:, 0] = polar_centered[:, 0] + (self._diameters[idx] / 2)
            self._full_xy = polar2cartesian(polar_centered) + center

        return self._full_xy

    @property
    def full_field_area(self):
        if self._full_area is None:
            self._full_area = spatial.ConvexHull(self.full_xy).volume

        return self._full_area
