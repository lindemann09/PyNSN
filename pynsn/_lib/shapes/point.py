__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from .abc_shape import ABCShape
from .picture_file import PictureFile
from ..coordinate import Coordinate


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
                PictureFile.check_attribute(attribute) is not None:
            raise NotImplementedError("Point arrays can not handle pictures.")

    def __repr__(self):
        return "Point(xy={}, attribute='{}')".format(self.xy, self.attribute)

    def distance(self, other):
        Coordinate.distance(self, other)

    @property
    def area(self):
        return 0

    @property
    def perimeter(self):
        return 0
