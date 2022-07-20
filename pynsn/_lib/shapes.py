__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math as _math
from .picture_file import PictureFile
from ..image._colour import Colour


class Point(object):

    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __repr__(self):
        return "Point(xy={})".format(self.xy)

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return Point(self.x * other, self.y * other)

    def __div__(self, other):
        return Point(self.x / other if other else self.x,
                     self.y / other if other else self.y)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        return self

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        return self

    def __imul__(self, other):
        self.x *= other
        self.y *= other
        return self

    def __idiv__(self, other):
        self.x /= other
        self.y /= other
        return self

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return self.x != other.x or self.y != other.y

    @property
    def xy(self):
        return self.x, self.y

    @xy.setter
    def xy(self, xy_tuple):
        self.x = xy_tuple[0]
        self.y = xy_tuple[1]

    @property
    def polar_radius(self):
        return _math.hypot(self.x, self.y)

    @polar_radius.setter
    def polar_radius(self, value):
        self.polar = (value, self.polar_angle)

    @property
    def polar_angle(self):
        return _math.atan2(self.y, self.x)

    @polar_angle.setter
    def polar_angle(self, value):
        self.polar = (self.polar_radius, value)

    @property
    def polar(self):
        """polar coordinate (radius, pos_angle) """
        return self.polar_radius, self.polar_angle

    @polar.setter
    def polar(self, rad_ang):
        """polar coordinate (radius, angle) """

        self.x = rad_ang[0] * _math.cos(rad_ang[1])
        self.y = rad_ang[0] * _math.sin(rad_ang[1])

    def distance(self, d):
        """Returns Euclidean distance to the another Coordinate. The function
        does not takes the size of an object into account.

        Parameters
        ----------
        d : Point

        Returns
        -------
        distance : float

        """

        return _math.hypot(self.x - d.x, self.y - d.y)


class ShapeAttribute(object):

    def __init__(self, attribute):
        self._attribute = None
        self.attribute = attribute  # setter

    @property
    def attribute(self):
        return self._attribute

    @attribute.setter
    def attribute(self, attr):
        """set attribute 
        
        Parameters
        ----------
        attr : anything
            can be, in principle, anything.
            If Colour or PictureFile, it will convert it to their string
            representation
        """
        if isinstance(attr, Colour):
            self._attribute = attr.colour
        elif isinstance(attr, PictureFile):
            self._attribute = attr.attribute
        else:
            self._attribute = attr

    def get_attribute_object(self):
        """Class instance of the attribute, if possible

        Returns
        -------
        rtn : attribute
            If attribute represents Colour or PictureFile it returns the instance
            of the respective class otherwise the previously defined  attribute
        """

        if isinstance(self._attribute, str):
            # check is color or picture
            col = Colour(self._attribute)
            if col.colour is not None:
                return col
            else:
                fl_name = PictureFile.check_attribute(self._attribute)
                if fl_name is not None:
                    return PictureFile(fl_name)
        return self._attribute


class Dot(Point, ShapeAttribute):

    def __init__(self, xy, diameter, attribute=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        xy : tuple of two numeric
        diameter : numeric
        attribute : attribute (string, optional)
        """
        if isinstance(attribute, PictureFile):
            raise NotImplementedError("Dot arrays can not handle pictures.")

        Point.__init__(self, x=xy[0], y=xy[1])
        ShapeAttribute.__init__(self, attribute)
        self.diameter = diameter

    def __repr__(self):
        return "Dot(xy={}, diameter={}, attribute='{}')".format(self.xy,
                                                                self.diameter, self.attribute)

    def distance(self, d):
        """Return Euclidean distance to the dot d. The function takes the
        diameter of the points into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float

        """

        return Point.distance(self, d) - ((self.diameter + d.diameter) / 2.0)

    @property
    def area(self):
        return _math.pi * (self.diameter ** 2) / 4.0

    @property
    def perimeter(self):
        return _math.pi * self.diameter


class Rectangle(Point, ShapeAttribute):

    def __init__(self, xy, size, attribute=None):
        """Initialize a Rectangle

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Rectangle can also consist of a picture

        Parameters
        ----------
        xy : tuple
            tuple of two numeric (default=(0, 0))
        size : tuple
            tuple of two numeric (default=(0, 0))
        attribute : attribute (string) or PictureFile
        """

        Point.__init__(self, x=xy[0], y=xy[1])
        ShapeAttribute.__init__(self, attribute)
        self.width, self.height = size

    def __repr__(self):
        return "Rectangle(xy={}, size={}, attribute='{}')".format(self.xy,
                                                                  self.size, self.get_attribute_str())

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

    def iter_edges(self, xy_tuple=False):
        """iterator over all four edges

        Returns
        -------
        iterator over Points or tuple (x,y)
        """
        for x in [Point(self.left, self.top),
                  Point(self.right, self.top),
                  Point(self.right, self.bottom),
                  Point(self.left, self.bottom)]:
            if xy_tuple:
                yield x.xy
            else:
                yield x

    @property
    def size(self):
        return self.width, self.height

    @size.setter
    def size(self, size):
        self.width, self.height = size

    def is_point_inside_rect(self, xy):
        return (self.left <= xy[0] <= self.right and
                self.top <= xy[1] <= self.bottom)

    def overlaps_with(self, rect):
        d = self.xy_distances(rect)
        return not (d[0] > 0 or d[1] > 0)

    @property
    def proportion(self):
        """Proportion of the rectangle (width/height)"""
        return self.width / self.height

    @property
    def area(self):
        return self.width * self.height

    @property
    def perimeter(self):
        return 2 * (self.width + self.height)

    @property
    def diagonal(self):
        return _math.sqrt(self.width ** 2 + self.height ** 2)

    def xy_distances(self, other):
        """return distances on both axes between rectangles. 0 indicates
        overlap off edges along that dimension.
        """
        assert isinstance(other, Rectangle)
        # overlaps in x or y
        pos_dist = abs(self.x - other.x), abs(self.y - other.y)
        max_overlap_dist = (self.width + other.width) / 2, (self.height + other.height) / 2
        if pos_dist[0] <= max_overlap_dist[0]:
            dx = 0
        else:
            dx = pos_dist[0] - max_overlap_dist[0]

        if pos_dist[1] <= max_overlap_dist[1]:
            dy = 0
        else:
            dy = pos_dist[1] - max_overlap_dist[1]

        return dx, dy

    def distance(self, other):
        """Return Euclidean distance to the dot d. The function takes the
        size of the rectangles into account.

        Parameters
        ----------
        other : Rectangle

        Returns
        -------
        distance : float

        """
        dx, dy = self.xy_distances(other=other)
        return _math.hypot(dx, dy)
