import numpy as np
from scipy import spatial
from .._lib.geometry import cartesian2polar, polar2cartesian
from .. import _arrays


class ConvexHullBaseClass(object):
    """convenient wrapper class for calculation of convex hulls
    """

    def __init__(self, xy):
        try:
            self._convex_hull = spatial.ConvexHull(xy)
        except spatial.QhullError:
            self._convex_hull = None

    @property
    def xy(self):
        if self._convex_hull is None:
            return np.array([])
        else:
            return self._convex_hull.points[self._convex_hull.vertices, :]

    @property
    def field_area(self):
        if self._convex_hull is None:
            return np.nan
        else:
            return self._convex_hull.volume

    @property
    def object_indices(self):
        return self._convex_hull.vertices

    def __eq__(self, other):
        """convex hulls are assumed to be identical if same field area is
        identical and convex comprises the same amount of points"""
        return len(self.xy) == len(other.xy) and \
               self.field_area == other.field_area


class ConvexHullPositions(ConvexHullBaseClass):
    """convenient wrapper class for calculation of convex hulls

    using just the positions (ignores the object size"
    """

    def __init__(self, object_array):
        assert isinstance(object_array, (_arrays.PointArray, _arrays.DotArray,
                                         _arrays.RectangleArray))
        ConvexHullBaseClass.__init__(self, xy=object_array.xy)


class ConvexHull(ConvexHullBaseClass):
    """convenient wrapper class for calculation of convex hulls

    using outer positions, able to handle rects
    """

    def __init__(self, object_array):
        assert isinstance(object_array, (_arrays.PointArray, _arrays.DotArray,
                                         _arrays.RectangleArray))

        if isinstance(object_array, _arrays.DotArray):
            # centered polar coordinates
            minmax = np.array((np.min(object_array.xy, axis=0),
                               np.max(object_array.xy, axis=0)))
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)  # center outer positions
            xy_centered = object_array.xy - center
            # outer positions
            polar_centered = cartesian2polar(xy_centered)
            polar_centered[:, 0] = polar_centered[:, 0] + (object_array.diameter / 2)
            outer_xy = polar2cartesian(polar_centered) + center

        elif isinstance(object_array, _arrays.RectangleArray):
            # simple solution: get all edges
            edges = []
            for rect in object_array.iter_objects():
                edges.extend([e.xy for e in rect.iter_edges()])
            outer_xy = np.array(edges)

        else:  # generic attribute array
            outer_xy = object_array.xy

        ConvexHullBaseClass.__init__(self, xy=outer_xy)
        self._is_rect_array = isinstance(object_array, _arrays.RectangleArray)

    @property
    def object_indices(self):
        if self._is_rect_array:
            return np.int16(self._convex_hull.vertices / 4)
        else:
            return self._convex_hull.vertices

    @property
    def object_indices_unique(self):
        return np.unique(self.object_indices)
