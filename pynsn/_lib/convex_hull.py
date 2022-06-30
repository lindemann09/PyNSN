import numpy as np
from scipy import spatial
from .geometry import cartesian2polar, polar2cartesian
from . import arrays

class ConvexHullBaseClass(object):
    """convenient wrapper class for calculation of convex hulls
    """

    def _initialize(self, xy):
        self._xy = xy
        try:
            self._convex_hull = spatial.ConvexHull(xy)
        except IndexError:
            self._convex_hull = None

    @property
    def xy(self):
        if self._convex_hull is None:
            return np.array([])
        else:
            return self._xy[self._convex_hull.vertices, :]

    @property
    def field_area(self):
        if self._convex_hull is None:
            return np.nan
        else:
            return self._convex_hull.volume

    @property
    def indices(self):
        return self._convex_hull.vertices


class ConvexHullPositions(ConvexHullBaseClass):
    """convenient wrapper class for calculation of convex hulls

    using just the positions (ignores the object size"
    """

    def __init__(self, object_array):
        arrays._check_generic_array(object_array)
        self._initialize(object_array.xy)


class ConvexHull(ConvexHullBaseClass):
    """convenient wrapper class for calculation of convex hulls

    using outer positions
    """

    def __init__(self, object_array):
        arrays._check_generic_array(object_array)

        if isinstance(object_array, arrays.DotArray):
            # centered polar coordinates
            minmax = np.array((np.min(object_array.xy, axis=0),
                               np.max(object_array.xy, axis=0)))
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)  # center outer positions
            xy_centered = object_array.xy - center
            # outer positions
            polar_centered = cartesian2polar(xy_centered)
            polar_centered[:, 0] = polar_centered[:, 0] + (object_array.diameters / 2)
            outer_xy = polar2cartesian(polar_centered) + center

        elif isinstance(object_array, arrays.RectangleArray):
            # simple solution: get all edges
            edges = []
            for rect in object_array.get():
                edges.extend([e.xy for e in rect.edges()])
            outer_xy = np.array(edges)

        else: # Generic object array
            outer_xy = object_array.xy

        self._initialize(outer_xy)
