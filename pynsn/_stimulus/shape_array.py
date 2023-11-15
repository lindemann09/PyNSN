"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"


from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .._shapes import Dot, Ellipse, Picture, Point2D, PolygonShape, Rectangle
from .._shapes import ellipse_geometry as ellipse_geo
from .._shapes.abc_shapes import AbstractCircularShape, AbstractShape
from ..errors import NoSolutionError
from ..random._rng import WalkAround
from .convex_hull import ConvexHull
from .target_area import TargetArea
from .._misc import delete_elements, IntOVector


class ShapeArray(object):
    """Numpy Optimizes Representation of Array of shapes"""

    def __init__(self) -> None:
        self._xy = np.empty((0, 2), dtype=np.float64)
        self._sizes = np.empty((0, 2), dtype=np.float64)
        self._shapes = []
        self._cache_polygons = None
        self._cache_ids = None
        self._convex_hull = None
        self._pos_changed = np.empty(0, dtype=int)
        self._size_changed = np.empty(0, dtype=int)

    def get_xy(self) -> NDArray[np.float_]:
        """array of positions"""
        return self._xy.copy()

    def set_xy(self, val: NDArray[np.float_]):
        # lasy update of shapes and polygons
        if val.shape != self._xy.shape:
            raise ValueError(
                f"xy has to be a numpy array with shape={self._xy.shape}")
        changes = np.any(self._xy != val, axis=1)
        if np.any(changes):
            self._pos_changed = np.unique(np.append(np.flatnonzero(changes),
                                                    self._pos_changed))
            self._xy = val
            self._clear_cached()

    def get_sizes(self) -> NDArray[np.float_]:
        """array of sizes"""
        return self._sizes.copy()

    def set_sizes(self, val: NDArray[np.float_]):
        # lasy update of shapes
        if val.shape != self._sizes.shape:
            raise ValueError(
                f"xy has to be a numpy array with shape={self._sizes.shape}")
        changes = np.any(self._sizes != val, axis=1)
        if np.any(changes):
            self._size_changed = np.unique(np.append(np.flatnonzero(changes),
                                                     self._size_changed))
            self._sizes = val
            self._clear_cached()

    def scale(self, factor: Union[float, np.float_]):
        """scale the size of all shapes"""
        if factor != 1:
            self.set_sizes(self._sizes * factor)

    @property
    def shapes(self) -> List[AbstractShape]:
        """list of the shapes"""
        if len(self._pos_changed) > 0 or len(self._size_changed) > 0:
            # update of shape list is needed, because size or pos changed
            for i in self._pos_changed:
                self._shapes[i].xy = self._xy[i, :]
            for i in self._size_changed:
                self._shapes[i].size = self._sizes[i, :]
            self._pos_changed = np.empty(0, dtype=int)
            self._size_changed = np.empty(0, dtype=int)

        return self._shapes

    @property
    def shape_types(self) -> NDArray[np.str_]:
        return np.array([s.name() for s in self._shapes])

    @property
    def polygons(self) -> NDArray[shapely.Polygon]:
        if self._cache_polygons is None:
            # build polygons
            poly = [s.polygon for s in self.shapes]
            self._cache_polygons = np.array(poly)
        return self._cache_polygons

    @property
    def attributes(self) -> List[Any]:
        return [s.attribute for s in self._shapes]

    @property
    def ids(self) -> Dict[str, NDArray]:
        """Dictionary of ids for the different shape types"""

        if self._cache_ids is None:
            self._cache_ids = {}
            for st in [Dot.name(), Rectangle.name(), Ellipse.name(), Picture.name(), PolygonShape.name()]:
                self._cache_ids[st] = np.flatnonzero(self.shape_types == st)
        return self._cache_ids

    def _clear_cached(self):
        """clears convex hull and ids"""
        self._convex_hull = None
        self._cache_ids = None
        self._cache_polygons = None

    def todict(self) -> dict:
        """Dict representation of the shape array
        """

        return {"shape_array": [x.todict() for x in self.shapes]}

    def shape_table_dict(self) -> dict:
        """Tabular representation of the array of the shapes.

        This representation can not deal with PolygonShapes. It"s useful to
        create Pandas dataframe or Arrow Tables.

        Examples
        --------
        >>> df_dict = stimulus.shape_table_dict()
        >>> df = pandas.DataFame(df_dict) # Pandas dataframe

        >>> table = pyarrow.Table.from_pydict(df_dict) # Arrow Table
        """

        if np.any(self.shape_types == PolygonShape.name()):
            raise RuntimeError("tabular shape representation can not deal with "
                               "PolygonShapes")
        d = {"type": self.shape_types.tolist(),
             "x": self._xy[:, 0].tolist(),
             "y": self._xy[:, 1].tolist(),
             "width": self._sizes[:, 0].tolist(),
             "height": self._sizes[:, 1].tolist(),
             "attributes": [str(x) for x in self.attributes]
             }
        return d

    def add(self, shape: AbstractShape):
        """add shape to the array"""
        if isinstance(shape, AbstractShape):
            self._shapes.append(shape)
            self._xy = np.append(self._xy, np.atleast_2d(shape.xy), axis=0)
            self._sizes = np.append(
                self._sizes, np.atleast_2d(shape.size), axis=0)
            self._clear_cached()
        else:
            raise TypeError(
                f"Can't add '{type(shape)}'. That's not a ShapeType."
            )

    def join_shapes(self, other: ShapeArray) -> None:
        """join with shapes of other array"""
        for x in other.shapes:
            self.add(x)

    def replace(self, index: int, shape: AbstractShape):
        self._shapes[index] = shape
        self._xy[index, :] = shape.xy
        self._sizes[index, :] = shape.size
        self._clear_cached()

    def delete(self, index: IntOVector) -> None:

        self._shapes = delete_elements(self._shapes, index)
        self._xy = np.delete(self._xy, index, axis=0)
        self._sizes = np.delete(self._sizes, index, axis=0)
        self._clear_cached()

    def clear(self):
        """ """
        self._shapes = []
        self._xy = np.empty((0, 2), dtype=np.float64)
        self._sizes = np.empty((0, 2), dtype=np.float64)
        self._clear_cached()

    def pop(self, index: Optional[int] = None) -> AbstractShape:
        """Remove and return item at index"""
        if index is None:
            index = len(self._shapes) - 1
        rtn = self.shapes[index]
        self.delete(index)
        return rtn

    def sort_by_excentricity(self):
        """Sort order fo the objects in the array by excentricity, that is, by the
        distance the center of the convex hull (from close to far)
        """
        d = self.distances(Point2D(self.convex_hull.centroid))
        idx = np.argsort(d)
        self._xy = self._xy[idx, :]
        self._sizes = self._sizes[idx]
        self._shapes = [self._shapes[i] for i in idx]
        self._clear_cached()

    def get_shapes(self, index: IntOVector) -> List[AbstractShape]:
        """Returns list with multiple selected objects

        Since shapes property has to be list (and not a numpy.array), this method
        allows slicing with a vector of ids (similar to numpy)

        Note
        ----
        Use `shapes`if you need to entire list of shapes.
        """
        if isinstance(index, int):
            ids = (index,)
        else:
            ids = index

        return [self.shapes[x] for x in ids]

    # def round_values(self, decimals: int = 0, int_type: type = np.int64) -> None: #FIXME rounding
    #     """rounds all values"""

    #     if decimals is None:
    #         return
    #     self._xy = round2(self._xy, decimals=decimals, int_type=int_type)
    #     self._dot_diameter = round2(self._dot_diameter, decimals=decimals,
    #                                 int_type=int_type)
    #     self._rect_sizes = round2(self._rect_sizes, decimals=decimals,
    #                               int_type=int_type)

    @property
    def n_objects(self) -> int:
        """number of shapes"""
        return len(self._shapes)

    @property
    def convex_hull(self) -> ConvexHull:
        """Convex hull polygon of the Shape Array"""
        if not isinstance(self._convex_hull, ConvexHull):
            self._convex_hull = ConvexHull(self.polygons)
        return self._convex_hull

    def csv(self, skip_columns: Optional[Sequence[str]] = None) -> str:
        """Comma-separated table representation the nsn stimulus

        Args:
            variable_names: if True first line include the variable names
            hash_column: if True hash will be included as first column
            skip_columns: list of column names that should not be exported
        Returns:
            CSV representation

        """

        d = self.todict()
        rtn = ",".join(d.keys()) + "\n"
        keys = list(d.keys())
        if skip_columns is not None:
            for sc in skip_columns:
                try:
                    keys.remove(sc)
                except ValueError:
                    pass

        for i in range(len(self._shapes)):
            row = ""
            for k in keys:
                row += f"{d[k][i]},"
            row = row[:-1] + "\n"  # replace comma
            rtn += row
        return rtn

    # def iter(self) -> Iterator:
    #     """iterate over all objects

    #     Parameters
    #     ----------
    #     indices: int or iterable of integer

    #     Notes
    #     -----
    #     To iterate all object you might all use the class iterator __iter__:
    #     >>> for obj in my_array:
    #     >>>    print(obj)
    #     """
    #     pass

    def distances(self, shape: Union[Point2D, AbstractShape]) -> NDArray[np.float_]:
        """distances of a shape or Point2D to the elements in the shape array"""

        if isinstance(shape, (Point2D, AbstractCircularShape)):
            # circular target shape
            rtn = np.full(len(self._shapes), np.nan)

            idx = self.ids[Dot.name()]
            if len(idx) > 0:
                # circular -> dots in shape array
                rtn[idx] = _distance_circ_dot_array(obj=shape,
                                                    dots_xy=self._xy[idx, :],
                                                    dots_diameter=self._sizes[idx, 0])
            idx = self.ids[Ellipse.name()]
            if len(idx) > 0:
                # circular -> ellipses in shape array
                rtn[idx] = _distance_circ_ellipse_array(obj=shape,
                                                        ellipses_xy=self._xy[idx, :],
                                                        ellipse_sizes=self._sizes[idx, :])
            # check if non-circular shapes are in shape_array
            idx = np.flatnonzero(np.isnan(rtn))
            if len(idx) > 0:
                if isinstance(shape, AbstractCircularShape):
                    rtn[idx] = shapely.distance(
                        shape.polygon, self.polygons[idx])
                else:
                    rtn[idx] = shapely.distance(
                        shape.xy_point, self.polygons[idx])
            return rtn

        else:
            # non-circular shape as target
            return shapely.distance(shape.polygon, self.polygons)

    def dwithin(self, shape: Union[Point2D, AbstractShape],  distance: float = 0) -> NDArray[np.bool_]:
        """Returns True for all elements of the array that are within the
        specified distance.

        Note
        -----
        Using this function is more efficient than computing the distance and comparing the result.
        """

        if isinstance(shape, (Point2D, AbstractCircularShape)):
            rtn = np.full(len(self._shapes), False)

            idx = self.ids[Dot.name()]
            if len(idx) > 0:
                # circular -> dots in shape array
                dists = _distance_circ_dot_array(obj=shape,
                                                 dots_xy=self._xy[idx, :],
                                                 dots_diameter=self._sizes[idx, 0])
                rtn[idx] = dists < distance

            idx = self.ids[Ellipse.name()]
            if len(idx) > 0:
                # circular -> ellipses in shape array
                dists = _distance_circ_ellipse_array(obj=shape,
                                                     ellipses_xy=self._xy[idx, :],
                                                     ellipse_sizes=self._sizes[idx, :])
                rtn[idx] = dists < distance

            # check if non-circular shapes are in shape_array
            idx = np.flatnonzero(np.isnan(rtn))
            if len(idx) > 0:
                if isinstance(shape, AbstractCircularShape):
                    rtn[idx] = shapely.dwithin(
                        shape.polygon, self.polygons[idx],  distance=distance)
                else:
                    rtn[idx] = shapely.dwithin(
                        shape.xy_point, self.polygons[idx],  distance=distance)

        else:
            # non-circular shape
            rtn = shapely.dwithin(
                shape.polygon, self.polygons, distance=distance)

        return rtn

    def contains_overlaps(self, min_distance: float = 0) -> bool:
        """Returns True for two or more elements overlap (i.e. taking
        into account the minimum distance).
        """
        for x in range(len(self._shapes)):
            if np.any(self.get_overlaps(x, min_distance)) > 0:
                return True
        return False

    def get_overlaps(self, index: int, min_distance: float = 0) -> NDArray[np.bool_]:
        """get overlaps with other objects. Ignores overlap with oneself."""
        overlaps = self.dwithin(self.shapes[index], distance=min_distance)
        overlaps[index] = False  # ignore overlap with oneself
        return overlaps

    @staticmethod
    def from_dict(the_dict: Dict[str, Any]) -> ShapeArray:
        """read shape array from dict"""
        ##
        raise NotImplementedError()  # FIXME

    def _random_free_position(self,
                              shape: AbstractShape,
                              min_distance: float,
                              target_area: TargetArea,
                              ignore_overlaps: bool,
                              max_iterations: int) -> AbstractShape:
        """returns the object at random free position

        raises exception if not found
        """

        cnt = 0
        while True:
            if cnt > max_iterations:
                raise NoSolutionError(
                    "Can't find a free position for this polygon")
            cnt += 1
            # propose a random position
            shape.xy = target_area.random_xy_inside_bounds()

            if not target_area.is_object_inside(shape):
                continue

            if ignore_overlaps:
                return shape
            else:
                # find overlaps
                overlaps = self.dwithin(shape, distance=min_distance)
                if not np.any(overlaps):
                    return shape

    def _fix_overlap(self,
                     index: int,
                     min_distance: float,
                     minimal_replacing: bool,
                     target_area: TargetArea,
                     max_iterations: int) -> int:
        """Move an selected object that overlaps to an free position in the
        neighbourhood.

        minimal_replacing: try to find a new random position is a neighbourhood,
            otherwise overlapping object will be randomly replaced anywere in the
            search area


        Returns
        -------
         0: if no overlaps exist
        -1: if object overlaps, but no new position could be found
         1: if object was replaced

        occupied space: see generator generate
        """

        if not np.any(self.get_overlaps(index, min_distance)):
            return 0  # no overlap

        target = self.shapes[index]

        if minimal_replacing:
            walk = WalkAround(target.xy)
            outside_cnt = 0
            while True:
                if walk.counter > max_iterations or outside_cnt > 20:
                    return -1  # can't find a free position

                target.xy = walk.next()
                if not target_area.is_object_inside(target):
                    outside_cnt += 1
                else:
                    outside_cnt = 0
                    overlaps = self.dwithin(target, distance=min_distance)
                    overlaps[index] = False  # ignore overlap with oneself
                    if not np.any(overlaps):
                        break  # place found
        else:
            # random position anywhere
            try:
                target = self._random_free_position(target,
                                                    min_distance=min_distance,
                                                    target_area=target_area,
                                                    ignore_overlaps=False,
                                                    max_iterations=max_iterations)
            except NoSolutionError:
                return -1

        self.replace(index, target)
        return 1


def _distance_circ_dot_array(obj: Union[Point2D, AbstractCircularShape],
                             dots_xy: NDArray[np.float_],
                             dots_diameter: NDArray[np.float_]) -> NDArray[np.float_]:
    """Distances circular shape or Point to multiple dots
    """
    d_xy = dots_xy - obj.xy
    if isinstance(obj, Point2D):
        circ_dia = 0
    elif isinstance(obj, Dot):
        circ_dia = obj.diameter
    elif isinstance(obj, Ellipse):
        circ_dia = ellipse_geo.diameter(
            size=np.atleast_2d(obj.size),
            theta=np.arctan2(d_xy[:, 1], d_xy[:, 0]))  # ellipse radius to each dot in the array
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(obj)}")

    # center dist - radius_a - radius_b
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - (dots_diameter + circ_dia) / 2


def _distance_circ_ellipse_array(obj: Union[Point2D, AbstractCircularShape],
                                 ellipses_xy: NDArray[np.float_],
                                 ellipse_sizes: NDArray[np.float_]) -> NDArray[np.float_]:
    """Distance circular shape or Point2D to multiple ellipses
    """
    d_xy = ellipses_xy - obj.xy
    theta = np.arctan2(d_xy[:, 1], d_xy[:, 0])
    # radii of ellipses in array to circ_shape
    ellipse_dia = ellipse_geo.diameter(size=ellipse_sizes, theta=theta)
    if isinstance(obj, Point2D):
        shape_dia = 0
    elif isinstance(obj, Dot):
        shape_dia = obj.diameter
    elif isinstance(obj, Ellipse):
        shape_dia = ellipse_geo.diameter(size=np.atleast_2d(obj.size),
                                         theta=theta)  # ellipse radius to each ellipse in the array
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(obj)}")

    # center dist - radius_a - radius_b
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - (ellipse_dia + shape_dia)/2
