from __future__ import annotations

from pynsn._lib import geometry

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Union, Any
import numpy as np
from numpy.typing import ArrayLike, NDArray

from .abc_shape import ABCShape
from .picture_file import PictureFile
from .._lib.coordinate import Coordinate
from .. import _shapes
from .._lib.spatial_relations import DotDot


class Dot(ABCShape):
    __slots__ = ("_diameter", )

    def __init__(self, xy: ArrayLike = (0, 0),
                 diameter: float = 0,
                 attribute: Any = None) -> None:
        """Initialize a dot

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        xy : tuple of two numeric (default=(0,0))
        diameter : numeric (default=(0,0))
        attribute : attribute (string, optional)
        """
        if isinstance(attribute, PictureFile) or \
                PictureFile.is_picture_attribute(attribute):
            raise NotImplementedError("Dot_arrays can not handle pictures.")

        super().__init__(xy, attribute)
        self._diameter = diameter

    def __repr__(self):
        return f"Dot(xy={self.xy}, diameter={self._diameter}, " \
            + "attribute = '{self.attribute}')"

    @property
    def diameter(self) -> float:
        return self._diameter

    @diameter.setter
    def diameter(self, value: float) -> None:
        self._diameter = value

    def distance(self, other: Union[_shapes.Dot, _shapes.Rectangle, Coordinate]) -> float:
        # inherited doc

        if isinstance(other, _shapes.Dot):
            rel = DotDot(a_xy=self._xy,  # dist-functions required 2D data
                         a_diameter=np.array([self._diameter]),
                         b_xy=other.xy,
                         b_diameter=np.array([other.diameter]))
            return rel.distances()[0]

        elif isinstance(other, _shapes.Rectangle):
            return other.distance(self)

        elif other.__class__ == Coordinate:
            rel = DotDot(a_xy=self._xy,
                         a_diameter=np.array([self._diameter]),
                         b_xy=other.xy,
                         b_diameter=np.zeros(1))
            return rel.distances()[0]

        raise NotImplementedError(f"distance to {type(other)} "
                                  + "is implemented.")

    @property
    def area(self) -> float:
        return np.pi * (self._diameter ** 2) / 4.0

    @property
    def perimeter(self) -> float:
        return np.pi * self._diameter

    def is_inside(self, other: Union[_shapes.Dot, _shapes.Rectangle]) -> bool:
        # FIXME remove this function

        if isinstance(other, _shapes.Dot):
            return Coordinate.distance(self, other) \
                < (other.diameter - self._diameter)/2.0

        elif isinstance(other, _shapes.Rectangle):
            left_buttom = self._xy - self._diameter / 2
            right_top = self._xy + self._diameter / 2
            ltrb = other.get_ltrb().flatten()
            if left_buttom[0] < ltrb[0] \
                    or right_top[1] > ltrb[1] \
                    or right_top[0] > ltrb[2] \
                    or left_buttom[1] < ltrb[3]:
                return False
            else:
                return True

    def rectangles_inside(self, xy: NDArray, sizes: NDArray) -> Union[np.bool_, NDArray[np.bool_]]:
        # inherited doc
        # FIXME not tested
        # [[left, top], [right, botton]]
        corners = geometry.corners(rect_xy=xy, rect_sizes=sizes)
        xy_dist = corners - np.atleast_3d(self.xy)  # type: ignore
        distances = np.hypot(xy_dist[:, 0, :], xy_dist[:, 1, :])
        comp = distances < self.diameter/2
        if comp.shape[0] > 1:
            return np.all(comp, axis=1)
        else:
            return np.all(comp)

    def dots_inside(self, xy: NDArray, diameters: NDArray) -> Union[np.bool_, NDArray[np.bool_]]:
        # inherited doc
        # FIXME not tested
        xy_diff = np.atleast_2d(xy) - np.atleast_2d(self.xy)
        comp = np.hypot(xy_diff[:, 0], xy_diff[:, 1]) < \
            (diameters - self.diameter)/2
        if comp.shape[0] > 1:
            return np.all(comp, axis=1)
        else:
            return np.all(comp)
