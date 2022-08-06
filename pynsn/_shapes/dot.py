from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math
from typing import Tuple, Union, Any

from .abc_shape import ABCShape
from .picture_file import PictureFile
from .._lib.coordinate import Coordinate
from .. import _shapes


class Dot(ABCShape):

    def __init__(self, xy: Tuple[float, float] =(0, 0),
                 diameter: float=1,
                 attribute: Any =None):
        """Initialize a dot

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        xy : tuple of two numeric
        diameter : numeric
        attribute : attribute (string, optional)
        """
        if isinstance(attribute, PictureFile) or \
                PictureFile.is_picture_attribute(attribute):
            raise NotImplementedError("Dot _arrays can not handle pictures.")

        super().__init__(xy, attribute)
        self.diameter = diameter

    def __repr__(self):
        return "Dot(xy={}, diameter={}, attribute='{}')".format(self.xy,
                                                                self.diameter, self.attribute)

    def distance(self, other: Union[_shapes.ShapeType, Coordinate]) -> float:
        # inherited doc

        if isinstance(other, _shapes.Dot):
            # distance between centers - radii
            return Coordinate.distance(self, other) \
                - ((self.diameter + other.diameter) / 2.0)

        elif isinstance(other, _shapes.Rectangle):
            return other.distance(self)

        elif isinstance(other, _shapes.Point) or other.__class__ == Coordinate:
            # distance between centers - radius
            return Coordinate.distance(self, other) \
                - (self.diameter / 2.0)

        raise NotImplementedError(f"distance to {type(other)} "
                                  + "is implemented.")

    @ property
    def area(self):
        return math.pi * (self.diameter ** 2) / 4.0

    @ property
    def perimeter(self):
        return math.pi * self.diameter

    def is_inside(self, other: Union[_shapes.ShapeType, Coordinate]) -> bool:
        # inherited doc
        if isinstance(other, _shapes.Dot):
            return Coordinate.distance(self, other) \
                < (other.diameter - self.diameter)/2.0

        elif isinstance(other, _shapes.Rectangle):
            r_self = self.diameter/2.0
            if self.x-r_self < other.left \
                    or self.x+r_self > other.right \
                    or self.y-r_self < other.bottom \
                    or self.y+r_self > other.top:
                return False
            else:
                return True

        elif isinstance(other, _shapes.Point) or other.__class__ == Coordinate:
            return False

        raise NotImplementedError("is_inside is not "
                                  f"implemented for {type(other)}.")
