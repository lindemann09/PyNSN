"""Shape Classes

Rectanglur shapes and Polygons are based on shapely.Polygon.
For performance reasons, circular shapes (Dot & Ellipse) are represented merely by
positions and radii. Polygons will be created if required only.

Note: Spatial module between circular shapes not on shapely.polygons.
"""

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from copy import deepcopy
from pathlib import Path
from typing import Any, Optional, Tuple, Union

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import Polygon

from ..types import Coord2DLike
from .abc_shapes import PointType, ShapeType, is_in_shape


class Rectangle(ShapeType):

    def __init__(self,
                 xy: Coord2DLike,
                 size: Coord2DLike,
                 attribute: Any = None
                 ):
        super().__init__(xy=xy, attribute=attribute)
        # make polygon
        l = xy[0] - size[0] / 2
        r = xy[0] + size[0] / 2
        t = xy[1] + size[1] / 2
        b = xy[1] - size[1] / 2
        self._polygon = Polygon(((l, b), (l, t), (r, t), (r, b)))
        shapely.prepare(self._polygon)  # FIXME needed?

    @property
    def polygon(self) -> Polygon:
        return self._polygon

    @property
    def size(self) -> NDArray:
        lb = shapely.get_coordinates(self._polygon)[0]  # left bottom
        rt = shapely.get_coordinates(self._polygon)[2]  # right top
        return rt - lb

    @property
    def proportion(self) -> float:
        """Proportion of the rectangle (width/height)"""
        return self.size[0] / self.size[1]

    @property
    def left_bottom(self) -> NDArray:
        """Returns (left, bottom) as ndarray (x,y)"""
        return shapely.get_coordinates(self._polygon)[0]

    @property
    def right_top(self) -> NDArray:
        """Returns (right, top) as ndarray (x,y)"""
        return shapely.get_coordinates(self._polygon)[2]

    @property
    def left_top(self) -> NDArray[np.float_]:
        """Returns (left, top) as ndarray (x,y)"""
        return shapely.get_coordinates(self._polygon)[1]

    @property
    def right_bottom(self) -> NDArray[np.float_]:
        """Returns (right, bottom) as ndarray (x,y)"""
        return shapely.get_coordinates(self._polygon)[3]

    @property
    def box(self) -> NDArray[np.float_]:
        """Returns (left, bottom, right, top) as NDArray (x0, y0, x1, y1)"""
        return np.append(
            shapely.get_coordinates(self._polygon)[0],
            shapely.get_coordinates(self._polygon)[2])

    def __repr__(self):
        return (f"Rectangle(xy={self._xy}, size={self.size}, "
                + f"attribute='{self._attribute}')")

    def todict(self) -> dict:
        d = super().todict()
        d.update({"size": self.size.tolist()})
        return d

    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> Rectangle:

        if copy_polygon:
            new = deepcopy(self)
            if new_xy is not None:
                new.xy = new_xy
            return new
        elif new_xy is None:
            return Rectangle(xy=self._xy, size=self.size, attribute=self._attribute)
        else:
            return Rectangle(xy=new_xy, size=self.size, attribute=self._attribute)

    def distance(self, shape: Union[PointType, ShapeType]) -> float:
        if isinstance(shape, PointType):
            return shapely.distance(self.polygon, shape.xy_point)
        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: Union[PointType, ShapeType], dist: float) -> bool:
        if isinstance(shape, PointType):
            return shapely.dwithin(self.polygon, shape.xy_point, distance=dist)
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=dist)

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        return is_in_shape(self, b=shape,
                           b_exterior_ring=shape_exterior_ring,
                           min_dist_boarder=min_dist_boarder)



class Picture(Rectangle):

    def __init__(self, xy: Coord2DLike, size: Coord2DLike,
                 path: Union[Path, str]) -> None:
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

        super().__init__(xy=xy, size=size, attribute=Path(path))

    def __repr__(self):
        return (f"Picture(xy={self._xy}, size={self.size}, " +
                f"path='{str(self.path)}')")

    @property
    def path(self) -> Path:
        return self._attribute

    def file_exists(self) -> bool:
        """Checks if the file exists"""
        return self.path.isfile()  # type: ignore

    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> Picture:

        if copy_polygon:
            new = deepcopy(self)
            if new_xy is not None:
                new.xy = new_xy
            return new
        elif new_xy is None:
            return Picture(xy=self._xy, size=self.size,
                           path=self.path)
        else:
            return Picture(xy=new_xy, size=self.size,
                           path=self.path)

    def todict(self) -> dict:
        return super().todict()

class PolygonShape(ShapeType):

    def __init__(self, polygon: Polygon, attribute: Any = None):
        ctr = polygon.centroid
        super().__init__(xy=(ctr.x, ctr.y), attribute=attribute)
        shapely.prepare(polygon)
        self._polygon = polygon

    @property
    def polygon(self) -> Polygon:
        return self._polygon

    @property
    def size(self) -> NDArray:
        b = shapely.bounds(self._polygon)  # l, b, r, t
        return b[2:4] - b[0:2]  # [bound width, bound height]

    def __repr__(self):
        return (f"PolygonShape(xy={self._xy}, size={self.size}, "
                + f"attribute='{self._attribute}')")

    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> PolygonShape:

        if not copy_polygon:
            raise RuntimeError(
                "copy_polygon = False not is possible for PolygonShape")

        new = deepcopy(self)
        if new_xy is not None:
            new.xy = new_xy
        return new

    def distance(self, shape: Union[PointType, ShapeType]) -> float:
        if isinstance(shape, PointType):
            return shapely.distance(self.polygon, shape.xy_point)
        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: Union[PointType, ShapeType], dist: float) -> bool:
        if isinstance(shape, PointType):
            return shapely.dwithin(self.polygon, shape.xy_point, distance=dist)
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=dist)

    def is_inside(self, shape: ShapeType,
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
