from numpy.typing import NDArray
import numpy as np
from shapely import MultiPolygon, Polygon, get_coordinates


class ConvexHull(object):
    """convenient wrapper class for calculating of convex hulls"""

    def __init__(self, array_of_polygons: NDArray[Polygon]) -> None:
        polys = MultiPolygon(array_of_polygons.tolist())
        self._ch_polygon = polys.convex_hull

    @property
    def polygon(self) -> Polygon:
        return self._ch_polygon

    @property
    def coordinates(self) -> NDArray[np.float_]:
        """Coordinates that shape the convex hull

        Returns:
            numpy array (x, y) of coordinates
        """
        return get_coordinates(self._ch_polygon)

    @property
    def area(self) -> float:
        """Area inside the convex hull"""
        return self._ch_polygon.area

    @property
    def centroid(self) -> NDArray[np.float_]:
        """Geometric Center of convex hull

        Returns:
            numpy array (x, y) of center position
        """
        return get_coordinates(self._ch_polygon.centroid)[0]
