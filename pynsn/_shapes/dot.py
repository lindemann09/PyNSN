from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math

from .abc_shape import ABCShape
from .picture_file import PictureFile
from .._lib.coordinate import Coordinate
from .. import _shapes


class Dot(ABCShape):

    def __init__(self, xy, diameter, attribute=None):
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

    def distance(self, other: _shapes.ShapeType) -> float:
        # inherited doc

        if isinstance(other, _shapes.Dot):
            return Coordinate.distance(self, other) \
                - ((self.diameter + other.diameter) / 2.0)

        elif isinstance(other, _shapes.Rectangle):
            return other.distance(self)

        elif isinstance(other, _shapes.Point):
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

    def is_inside(self, other: _shapes.ShapeType) -> bool:
        # inherited doc
        if isinstance(other, _shapes.Dot):
            pass

        elif isinstance(other, _shapes.Rectangle):
            pass

        elif isinstance(other, _shapes.Point):
            return False

        raise NotImplementedError("is_inside is not "
                                  "implemented for {}.".format(type(other)))
