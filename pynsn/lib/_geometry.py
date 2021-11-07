
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math
import numpy as np



def polar2cartesian(polar):
    """polar is an 2d-array representing polar coordinates (radius, angle)"""
    polar = np.array(polar)
    return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                     polar[:, 0] * np.sin(polar[:, 1])]).T

def cartesian2polar(xy, radii_only=False):
    """polar coordinates (radius, angle)

    if only radii required you may consider radii_only=True for faster
    processing
    """
    xy = np.array(xy)
    radii = np.hypot(xy[:, 0], xy[:, 1])
    if radii_only:
        return radii
    else:
        return np.array([radii, np.arctan2(xy[:, 1], xy[:, 0])]).T


def lines_intersect(line1, line2):
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
def distance_between_edge_and_point(edge, point):
    # edge is a tuple of 2d coordinates
    if point_faces_edge(edge, point):
        area=triangle_area_at_points(edge[0], edge[1], point)
        base=edge[0].distance(edge[1])
        height=area/(0.5*base)
        return height
    return min(point.distance(edge[0]), point.distance(edge[1]))

def triangle_area_at_points(p1, p2, p3):
    # p1, p2, p3: 2d Coordinates
    a=p1.distance(p2)
    b=p2.distance(p3)
    c=p1.distance(p3)
    s=(a+b+c)/float(2)
    area=math.sqrt(s*(s-a)*(s-b)*(s-c))
    return area

# Finds angle using cos law
def angle(a, b, c):
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
def point_faces_edge(edge, point):
    a=edge[0].distance(edge[1])
    b=edge[0].distance(point)
    c=edge[1].distance(point)
    ang1, ang2 = angle(b, a, c), angle(c, a, b)
    if ang1>math.pi/2 or ang2>math.pi/2:
        return False
    return True
