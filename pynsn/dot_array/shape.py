
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math
from ..lib.coordinate2D import Coordinate2D
from .item_attributes import ItemAttributes


class Dot(Coordinate2D):  # TODO becomes maybe an item

    def __init__(self, x=0, y=0, diameter=1, attributes=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        diameter : numeric (default=1)
        attributes : ItemAttributes
        """

        Coordinate2D.__init__(self, x=x, y=y)
        self.diameter = diameter
        if attributes is None:
            self.attributes = ItemAttributes()
        elif not isinstance(attributes, ItemAttributes):
            raise TypeError("features must be a ItemFeatures, not {}".format(type(attributes).__name__))
        else:
            self.attributes = attributes

    def distance(self, d):
        """Return Euclidean distance to the dot d. The function takes the
        diameter of the dots into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float

        """

        return Coordinate2D.distance(self, d) - \
               ((self.diameter + d.diameter) / 2.0)

    @property
    def area(self):
        return math.pi * (self.diameter ** 2) / 4.0

    @property
    def perimeter(self):
        return math.pi * self.diameter


class Rectangle(Coordinate2D):

    def __init__(self, center_x=0, center_y=0, width=0, height=0,
                 attributes=None):
        """Initialize a Rectangle

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        width : numeric (default=1)
        height : numeric (default=1)
        """

        Coordinate2D.__init__(self, x=center_x, y=center_y)
        if attributes is None:
            self.attributes = ItemAttributes()
        else:
            self.attributes = attributes

        self.height = height
        self.width = width

    @property
    def left(self):
        return self.x - 0.5 * self.width

    @property
    def top(self):
        return self.y + 0.5 * self.height

    @property
    def right(self):
        return self.x + 0.5 * self.width

    @property
    def bottom(self):
        return self.y - 0.5 * self.height

    def is_point_inside_rect(self, xy):
        return (self.left <= xy[0] <= self.right and
                self.top <= xy[1] <= self.bottom)

    def overlaps_with(self, rect):
        d = self.xy_distances(rect)
        return not(d[0]>0 or d[1]>0)

    @property
    def area(self):
        return self.width * self.height

    @property
    def perimeter(self):
        return 2 * (self.width + self.height)

    def xy_distances(self, other):
        """return distances on both axes between rectangles. 0 indicates
        overlap off edges along that dimension.
        """

        assert isinstance(other, Rectangle)
        #overlaps in x or y
        center_dist = abs(self.x - other.x), abs(self.y - other.y)
        max_overlap_dist = (self.width + other.width)/2, (self.height + other.height)/2
        if  center_dist[0] <= max_overlap_dist[0]:
            dx = 0
        else:
            dx = center_dist[0] - max_overlap_dist[0]

        if  center_dist[1] <= max_overlap_dist[1]:
            dy = 0
        else:
            dy = center_dist[1] - max_overlap_dist[1]

        return dx, dy


    def distance(self, other):
        """Return Euclidean distance to the dot d. The function takes the
        size of the rectangles into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float

        """

        dx, dy = self.xy_distances(other=other)
        return math.hypot(dx, dy)
