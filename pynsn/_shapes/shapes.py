"""Shape Classes

Rectangular shapes and Polygons are based on shapely.Polygon.
For performance reasons, circular shapes (Dot & Ellipse) are represented merely by
positions and radii. Polygons will be created if required only.

Note: Spatial module between circular shapes not on shapely.polygons.
"""

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from pathlib import Path
from typing import Optional, Union

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import Polygon

from .abc_shapes import (AbstractPoint, AbstractShape, Coord2DLike,
                         is_in_shape, AttributeType)


class Rectangle(AbstractShape):

    def __init__(self,
                 size: Coord2DLike,
                 xy: Coord2DLike = (0, 0),
                 attribute: AttributeType = None
                 ):
        """A Rectangle Shape"""
        super().__init__(size=size, xy=xy, attribute=attribute)

    @property
    def polygon(self) -> Polygon:
        if self._polygon is None:
            # make polygon
            l = self._xy[0] - self._size[0] / 2
            r = self._xy[0] + self._size[0] / 2
            t = self._xy[1] + self._size[1] / 2
            b = self._xy[1] - self._size[1] / 2
            self._polygon = Polygon(((l, b), (l, t), (r, t), (r, b)))
            shapely.prepare(self._polygon)

        return self._polygon

    @property
    def proportion(self) -> float:
        """Proportion of the rectangle (width/height)"""
        return self.size[0] / self.size[1]

    @property
    def left_bottom(self) -> NDArray:
        """Returns (left, bottom) as ndarray (x,y)"""
        return shapely.get_coordinates(self.polygon)[0]

    @property
    def right_top(self) -> NDArray:
        """Returns (right, top) as ndarray (x,y)"""
        return shapely.get_coordinates(self.polygon)[2]

    @property
    def left_top(self) -> NDArray[np.float_]:
        """Returns (left, top) as ndarray (x,y)"""
        return shapely.get_coordinates(self.polygon)[1]

    @property
    def right_bottom(self) -> NDArray[np.float_]:
        """Returns (right, bottom) as ndarray (x,y)"""
        return shapely.get_coordinates(self.polygon)[3]

    @property
    def box(self) -> NDArray[np.float_]:
        """Returns (left, bottom, right, top) as NDArray (x0, y0, x1, y1)"""
        return np.append(
            shapely.get_coordinates(self.polygon)[0],
            shapely.get_coordinates(self.polygon)[2])

    def __repr__(self):
        return (f"Rectangle(xy={self._xy}, size={self.size}, "
                + f"attribute='{self._attribute}')")

    def todict(self) -> dict:
        d = super().todict()
        d.update({"size": self.size.tolist()})
        return d

    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        if isinstance(shape, AbstractPoint):
            return shapely.distance(self.polygon, shape.xy_point)
        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: Union[AbstractPoint, AbstractShape], distance: float) -> bool:
        if isinstance(shape, AbstractPoint):
            return shapely.dwithin(self.polygon, shape.xy_point, distance=distance)
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=distance)

    def is_inside(self, shape: AbstractShape,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        return is_in_shape(self, b=shape,
                           b_exterior_ring=shape_exterior_ring,
                           min_dist_boarder=min_dist_boarder)


class Picture(Rectangle):

    def __init__(self,
                 size: Coord2DLike,
                 path: Union[Path, str],
                 xy: Coord2DLike = (0, 0)) -> None:
        """Initialize a Picture

        Rectangle can also consist of a picture

        Parameters
        ----------
        xy : tuple
            tuple of two numeric
        size : tuple
            tuple of two numeric
        path : pathlib.path or str
            path the picture file
        """

        super().__init__(size=size, xy=xy, attribute=Path(path))

    def __repr__(self):
        return (f"Picture(xy={self._xy}, size={self.size}, " +
                f"path='{str(self.path)}')")

    @property
    def path(self) -> Path:
        return Path(self._attribute)  # type: ignore

    def file_exists(self) -> bool:
        """Checks if the file exists"""
        return self.path.is_file()


class PolygonShape(AbstractShape):

    def __init__(self, polygon: Polygon, attribute: AttributeType = None):
        if not isinstance(polygon, Polygon):
            raise TypeError(
                f"Polygon has to be a shapely.Polygon and not a {type(polygon)}")

        ctr = polygon.centroid
        b = shapely.bounds(polygon)  # l, b, r, t
        size = b[2:4] - b[0:2]  # [bound width, bound height]
        shapely.prepare(polygon)

        super().__init__(size=size, xy=(ctr.x, ctr.y), attribute=attribute)
        self._polygon = polygon

    @property
    def polygon(self) -> Polygon:
        return self._polygon

    def __repr__(self):
        return (f"PolygonShape(xy={self._xy}, size={self.size}, "
                + f"attribute='{self._attribute}')")

    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        if isinstance(shape, AbstractPoint):
            return shapely.distance(self.polygon, shape.xy_point)
        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: Union[AbstractPoint, AbstractShape], distance: float) -> bool:
        if isinstance(shape, AbstractPoint):
            return shapely.dwithin(self.polygon, shape.xy_point, distance=distance)
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=distance)

    def is_inside(self, shape: AbstractShape,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        return is_in_shape(self, b=shape,
                           b_exterior_ring=shape_exterior_ring,
                           min_dist_boarder=min_dist_boarder)

    def todict(self) -> dict:
        d = super().todict()
        del d["xy"]
        d.update({"wkt": shapely.to_wkt(self.polygon)})
        return d
