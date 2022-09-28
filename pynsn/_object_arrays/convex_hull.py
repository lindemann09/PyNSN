from __future__ import annotations

from abc import ABCMeta, abstractmethod

import numpy as np
from numpy.typing import NDArray
from scipy import spatial

from .._lib.geometry import cartesian2polar, polar2cartesian, center_of_coordinates
from . import ObjectArrayType
from .dot_array import DotArray
from .rectangle_array import RectangleArray


class ABCConvexHull(metaclass=ABCMeta):
    """convenient wrapper class for calculating of convex hulls

    using just the positions (ignores the object size"
    """

    @abstractmethod
    def __init__(self) -> None:
        self._convex_hull = None

    def _try_convex_hull(self, xy):
        # pylint: disable=E1101
        try:
            self._convex_hull = spatial.ConvexHull(xy)
        except spatial.QhullError:
            self._convex_hull = None

    @property
    def xy(self) -> NDArray:
        '''Coordinates that shape the convex hull

        Returns:
            np.array of coordinates
        '''
        if self._convex_hull is None:
            return np.array([])
        else:
            return self._convex_hull.points[self._convex_hull.vertices, :]

    @property
    def field_area(self) -> float:
        '''Area inside the convex hull
        '''
        if self._convex_hull is None:
            return 0
        else:
            return self._convex_hull.volume

    @property
    def _point_indices(self) -> NDArray:
        '''Indices of the points (xy) that shape the convex hull

        Returns:
            np.array with indices
        '''
        if self._convex_hull is None:
            return np.array([])
        else:
            return self._convex_hull.vertices

    def __eq__(self, other):
        """convex hulls are assumed to be identical if same field area is
        identical and convex comprises the same amount of points"""
        return len(self.xy) == len(other.xy) and \
            self.field_area == other.field_area

    @property
    def center(self) -> NDArray:
        """Center of convex hull

        Returns:
            Coordinate of center position
        """

        return center_of_coordinates(self.xy)


class ConvexHullPositions(ABCConvexHull):

    def __init__(self, object_array: ObjectArrayType) -> None:
        super().__init__()
        assert isinstance(object_array, (RectangleArray, DotArray))
        self._try_convex_hull(xy=object_array.xy)


class ConvexHull(ABCConvexHull):
    """convenient wrapper class for calculation of convex hulls

    using outer positions, able to handle rects
    """

    def __init__(self, object_array: ObjectArrayType) -> None:

        super().__init__()
        assert isinstance(object_array, (RectangleArray, DotArray))
        self._is_rect_array = isinstance(object_array, RectangleArray)
        if self._is_rect_array:
            if len(object_array.xy) < 2:
                return  # no convex hull possible
        else:
            # object point array at least 3 points
            if len(object_array.xy) < 2:
                return  # no convex hull possible

        if isinstance(object_array, DotArray):
            if len(object_array.xy) < 2:
                return  # no convex hull possible
            # centered polar coordinates
            minmax = np.array((np.min(object_array.xy, axis=0),
                               np.max(object_array.xy, axis=0)))
            # center outer positions
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)
            xy_centered = object_array.xy - center
            # outer positions
            polar_centered = cartesian2polar(xy_centered)
            polar_centered[:, 0] = polar_centered[:, 0] + \
                (object_array.diameter / 2)
            outer_xy = polar2cartesian(polar_centered) + center

        elif isinstance(object_array, RectangleArray):
            # simple solution: get all corners
            corners = []
            for rect in object_array.iter():
                corners.extend([e.xy for e in rect.iter_corners()])
            outer_xy = np.array(corners)

        else:  # generic attribute array
            outer_xy = object_array.xy

        self._try_convex_hull(xy=outer_xy)

    @property
    def object_indices(self):
        '''Indices of the objects of the associated array are part of the
        convex hull

        Returns:
            np.array with the object ids
        '''
        if self._is_rect_array:
            rtn = np.int16(self._point_indices / 4)
            return np.unique(rtn)
        else:
            return self._point_indices