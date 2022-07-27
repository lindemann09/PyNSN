__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math

from .abc_shape import ABCShape
from .picture_file import PictureFile
from .._lib.coordinate import Coordinate


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

    def distance(self, other):
        """Return Euclidean distance to the dot d. The function takes the
        diameter of the points into account.

        Parameters
        ----------
        other : Dot

        Returns
        -------
        distance : float

        """

        return Coordinate.distance(self, other) - \
               ((self.diameter + other.diameter) / 2.0)

    @property
    def area(self):
        return math.pi * (self.diameter ** 2) / 4.0

    @property
    def perimeter(self):
        return math.pi * self.diameter


