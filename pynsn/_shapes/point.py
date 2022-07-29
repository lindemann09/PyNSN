from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from .abc_shape import ABCShape
from .picture_file import PictureFile
from .._lib.coordinate import Coordinate
from .. import _shapes


class Point(ABCShape):

    def __init__(self, xy, attribute=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        xy : tuple of two numeric
        attribute : attribute (string, optional)
        """

        super().__init__(xy, attribute)
        if isinstance(attribute, PictureFile) or \
                PictureFile.is_picture_attribute(attribute):
            raise NotImplementedError("Point _arrays can not handle pictures.")

    def __repr__(self):
        return "Point(xy={}, attribute='{}')".format(self.xy, self.attribute)

    @property
    def area(self):
        return 0

    @property
    def perimeter(self):
        return 0

    def distance(self, other: _shapes.ShapeType) -> float:
        # inherited doc

        if isinstance(other, _shapes.Point):
            return Coordinate.distance(self, other)
        elif isinstance(other, (_shapes.Dot, _shapes.Rectangle)):
            return other.distance(self)

        raise NotImplementedError(f"distance to {type(other)} "
                                  + "is implemented.")

    def is_inside(self, other: _shapes.ShapeType) -> bool:
        # inherited doc
        if isinstance(other, _shapes.Dot):
            pass

        elif isinstance(other, _shapes.Rectangle):
            pass

        elif isinstance(other, _shapes.Point):
            pass

        raise NotImplementedError("is_inside is not "
                                  "implemented for {}.".format(type(other)))
