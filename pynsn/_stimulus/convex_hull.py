import numpy as np
from numpy.typing import NDArray
from shapely import MultiPolygon, Polygon, get_coordinates

from .._lib.geometry import Coord2D


class ConvexHull(object):
    """convenient wrapper class for calculating of convex hulls"""

    def __init__(self, array_of_polygons: NDArray[Polygon]) -> None:
        polys = MultiPolygon(list(array_of_polygons))
        self._ch_polygon = polys.convex_hull

    @property
    def polygon(self) -> Polygon:
        return self._ch_polygon

    @property
    def coordinates(self) -> NDArray:
        """Coordinates that shape the convex hull

        Returns:
            np.array of coordinates
        """
        return get_coordinates(self._ch_polygon)

    @property
    def area(self) -> float:
        """Area inside the convex hull"""
        return self._ch_polygon.area

    @property
    def centroid(self) -> Coord2D:
        """Geometric Center of convex hull

        Returns:
            Coordinate of center position
        """
        cent = self._ch_polygon.centroid
        return cent.x, cent.y
