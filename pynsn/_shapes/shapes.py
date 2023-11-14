"""Shape Classes

Rectangular shapes and Polygons are based on shapely.Polygon.
For performance reasons, circular shapes (Dot & Ellipse) are represented merely by
positions and radii. Polygons will be created if required only.

Note: Spatial module between circular shapes not on shapely.polygons.
"""

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from copy import deepcopy
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import Polygon, affinity

from .abc_shapes import AbstractPoint, AbstractShape, Coord2DLike, is_in_shape


class Rectangle(AbstractShape):

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

    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        if isinstance(shape, AbstractPoint):
            return shapely.distance(self.polygon, shape.xy_point)
        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: Union[AbstractPoint, AbstractShape], dist: float) -> bool:
        if isinstance(shape, AbstractPoint):
            return shapely.dwithin(self.polygon, shape.xy_point, distance=dist)
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=dist)

    def is_inside(self, shape: AbstractShape,
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
        return self._attribute # type: ignore

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


class PolygonShape(AbstractShape):

    def __init__(self, polygon: Polygon, attribute: Any = None):
        if not isinstance(polygon, Polygon):
            raise TypeError(f"Polygon has to be a shapely.Polygon and not a {type(polygon)}")

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

    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        if isinstance(shape, AbstractPoint):
            return shapely.distance(self.polygon, shape.xy_point)
        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: Union[AbstractPoint, AbstractShape], dist: float) -> bool:
        if isinstance(shape, AbstractPoint):
            return shapely.dwithin(self.polygon, shape.xy_point, distance=dist)
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=dist)

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

    def copy_scaled(self,
                    new_size:Coord2DLike,
                    new_xy: Optional[Coord2DLike] = None):
        """returns a copy of the polygon scaled to a different size and optionally
        new position"""
        new_size = np.asarray(new_size)
        if len(new_size) != 2:
            raise ValueError("new_size has be an list of two numerals (width, height)")
        fact = new_size / self.size
        rtn = PolygonShape(
            affinity.scale(self._polygon, xfact=fact[0], yfact=fact[1]))
        if new_xy is not None:
            rtn.xy = new_xy
        return rtn
