from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Union, Any
import numpy as np
from numpy.typing import ArrayLike

from .abc_shape import ABCShape
from .picture_file import PictureFile
from .._lib.coordinate import Coordinate
from .. import _shapes
from .._lib.geometry import dist_dots
from .._lib.misc import numpy_array_2d


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
            rtn = dist_dots(a_xy=numpy_array_2d(self.xy),  # dist-functions required 2D data
                            a_diameter=self._diameter,
                            b_xy=other.xy, b_diameter=other.diameter)
            return rtn[0]

        elif isinstance(other, _shapes.Rectangle):
            return other.distance(self)

        elif other.__class__ == Coordinate:
            rtn = dist_dots(a_xy=numpy_array_2d(self.xy),
                            a_diameter=self._diameter,
                            b_xy=other.xy,
                            b_diameter=0)
            return rtn[0]

        raise NotImplementedError(f"distance to {type(other)} "
                                  + "is implemented.")

    @ property
    def area(self) -> float:
        return np.pi * (self._diameter ** 2) / 4.0

    @ property
    def perimeter(self) -> float:
        return np.pi * self._diameter

    def is_inside(self, other: Union[_shapes.Dot, _shapes.Rectangle, Coordinate]) -> bool:
        # inherited doc
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

        elif other.__class__ == Coordinate:
            return False

        raise NotImplementedError("is_inside is not "
                                  f"implemented for {type(other)}.")
