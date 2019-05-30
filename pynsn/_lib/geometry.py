"""
Dot Array
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math
from .item_attributes import ItemAttributes, ItemAttributesList

class Coordinate2D(object):

    def __init__(self, x=0, y=0):
        self._x = x
        self._y = y
        self._polar_radius = None
        self._polar_angle = None

    @property
    def x(self):
        if self._x is None:
            self._x = self._polar_radius * math.cos(self._polar_angle)
        return self._x

    @x.setter
    def x(self, value):
        self.xy = (value, self.y)

    @property
    def y(self):
        if self._y is None:
            self._y = self._polar_radius * math.sin(self._polar_angle)
        return self._y

    @y.setter
    def y(self, value):
        self.xy = (self.x, value)

    @property
    def xy(self):
        return (self.x, self.y)

    @xy.setter
    def xy(self, xy_tuple):
        self._x = xy_tuple[0]
        self._y = xy_tuple[1]
        self._polar_radius = None
        self._polar_angle = None

    @property
    def polar_radius(self):
        if self._polar_radius is None:
            self._polar_radius = math.hypot(self._x, self._y)
        return self._polar_radius

    @polar_radius.setter
    def polar_radius(self, value):
        self.polar = (value, self.polar_angle)

    @property
    def polar_angle(self):
        if self._polar_angle is None:
            self._polar_angle = math.atan2(self._y, self._x)
        return self._polar_angle

    @polar_angle.setter
    def polar_angle(self, value):
        self.polar = (self.polar_radius, value)

    @property
    def polar(self):
        """polar coordinate (radius, pos_angle) """
        return (self.polar_radius, self.polar_angle)

    @polar.setter
    def polar(self, rad_ang):
        """polar coordinate (radius, angle) """

        self._polar_radius = rad_ang[0]
        self._polar_angle = rad_ang[1]
        self._x = None
        self._y = None

    def distance(self, d):
        """Return Euclidean distance to the another Coordinate. The function take the
        diameter of the points into account.

        Parameters
        ----------
        d : Coordinate2D

        Returns
        -------
        distance : float

        """

        return math.hypot(self.x - d.x, self.y - d.y)


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
            self.attributes = ItemAttributes(colour=None, picture=None)
        elif not isinstance(attributes, ItemAttributes):
            raise TypeError("features must be a ItemFeatures, not {}".format(type(attributes).__name__))
        else:
            self.attributes = attributes

    def distance(self, d):
        """Return Euclidean distance to the dot d. The function take the
        diameter of the points into account.

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

    def __init__(self, center_x=0, center_y=0, width=0, height=0, features=None):
        """Initialize a point

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
        if features is None:
            self.features = ItemAttributesList(colour=None)
        else:
            self.features = features

        self.height = height
        self.width = width

    @property
    def left(self):
        return self._x - 0.5 * self.width

    @property
    def top(self):
        return self._y + 0.5 * self.height

    @property
    def right(self):
        return self._x + 0.5 * self.width

    @property
    def bottom(self):
        return self._y - 0.5 * self.height


    def iter_edges(self):
        yield self.left, self.top
        yield self.right, self.top
        yield self.right, self.bottom
        yield self.left, self.bottom

    def is_point_inside_rect(self, xy):
        return (self.left <= xy[0] <= self.right and
                self.top <= xy[1] <= self.bottom)

    def overlaps_with(self, rect):
        for corner in rect.iter_edges():
            if self.is_point_inside_rect(corner):
                return True
        for corner in self.iter_edges():
            if rect.is_point_inside_rect(corner):
                return True
        return False

    def distance(self, rect):
        """Return Euclidean distance to other rect. The function take the
        diameter of the points into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float or -1 if overlapping

        """

        # 1. see if they overlap
        if self.overlaps_with(rect):
            return 0

        # 2. draw a line between rectangles
        line = (self.xy, rect.xy)

        # 3. find the two edges that intersect the line
        edge1 = None
        edge2 = None
        for edge in self.iter_edges():
            if _lines_intersect(edge, line):
                edge1 = edge
                break
        for edge in rect.iter_edges():
            if _lines_intersect(edge, line):
                edge2 = edge
                break
        assert edge1
        assert edge2

        # 4. find shortest distance between these two edges
        distances = [
            _distance_between_edge_and_point(edge1, edge2[0]),
            _distance_between_edge_and_point(edge1, edge2[1]),
            _distance_between_edge_and_point(edge2, edge1[0]),
            _distance_between_edge_and_point(edge2, edge1[1]),
        ]

        return min(distances)

    @property
    def area(self):
        return self.width * self.height

    @property
    def perimeter(self):
        return 2 * (self.width + self.height)


def _lines_intersect(line1, line2):
    # lines_overlap_on_x_axis
    x1, x2 = line1[0].x, line1[1].x
    x3, x4 = line2[0].x, line2[1].x
    e1_left, e1_right = min(x1, x2), max(x1, x2)
    e2_left, e2_right = min(x3, x4), max(x3, x4)
    lines_overlap_on_x_axis = (e1_left >= e2_left and e1_left <= e2_right) or \
            (e1_right >= e2_left and e1_right <= e2_right) or \
            (e2_left >= e1_left and e2_left <= e1_right) or \
            (e2_right >= e1_left and e2_right <= e1_right)

    # _lines_overlap_on_y_axis
    y1, y2 = line1[0].y, line1[1].y
    y3, y4 = line2[0].y, line2[1].y
    e1_top, e1_bot = min(y1, y2), max(y1, y2)
    e2_top, e2_bot = min(y3, y4), max(y3, y4)
    lines_overlap_on_y_axis = (e1_top >= e2_top and e1_top <= e2_bot) or \
           (e1_bot >= e2_top and e1_bot <= e2_bot) or \
           (e2_top >= e1_top and e2_top <= e1_bot) or \
           (e2_bot >= e1_top and e2_bot <= e1_bot)

    return lines_overlap_on_x_axis and lines_overlap_on_y_axis


# Gives distance if the point is facing edge, else False
def _distance_between_edge_and_point(edge, point):
    # edge is a tuple of 2d coordinates
    if _point_faces_edge(edge, point):
        area=_triangle_area_at_points(edge[0], edge[1], point)
        base=edge[0].distance(edge[1])
        height=area/(0.5*base)
        return height
    return min(point.distance(edge[0]), point.distance(edge[1]))

def _triangle_area_at_points(p1, p2, p3):
    # p1, p2, p3: 2d Coordinates
    a=p1.distance(p2)
    b=p2.distance(p3)
    c=p1.distance(p3)
    s=(a+b+c)/float(2)
    area=math.sqrt(s*(s-a)*(s-b)*(s-c))
    return area

# Finds angle using cos law
def _angle(a, b, c):
    divid=float(a**2+b**2-c**2)
    divis=(2*a*b)
    if (divis)>0:
        result=divid/divis
        if result<=1.0 and result>=-1.0:
            return math.acos(result)
        return 0
    else:
        return 0

# Checks if point faces edge
def _point_faces_edge(edge, point):
    a=edge[0].distance(edge[1])
    b=edge[0].distance(point)
    c=edge[1].distance(point)
    ang1, ang2 = _angle(b, a, c), _angle(c, a, b)
    if ang1>math.pi/2 or ang2>math.pi/2:
        return False
    return True
