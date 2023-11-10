"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"


from collections import OrderedDict
from copy import deepcopy
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .._shapes import (CircularShapeType, Dot, Ellipse, Picture, Point2D,
                       PolygonShape, Rectangle, ShapeType)
from .._shapes import ellipse_geometry as ellipse_geo
from ..random._rng import WalkAround
from ..types import IntOVector, NoSolutionError
from .convex_hull import ConvexHull
from .target_area import TargetArea


class ShapeArray(object):
    """Numpy Optimizes Representation of Array of shapes"""

    def __init__(self) -> None:
        self._xy = np.empty((0, 2), dtype=np.float64)
        self._sizes = np.empty((0, 2), dtype=np.float64)
        self._polygons = np.empty(0, dtype=object)
        self._attributes = np.empty(0, dtype=object)
        self._types = np.empty(0, dtype=str)
        self._ids = {}
        self._convex_hull = None
        self._update_ids()

    @property
    def xy(self) -> NDArray:
        return self._xy

    @property
    def sizes(self) -> NDArray:
        return self._sizes

    @property
    def shape_types(self) -> NDArray[np.str_]:
        return self._types

    @property
    def polygons(self) -> NDArray[shapely.Polygon]:
        return self._polygons

    @property
    def attributes(self) -> NDArray:
        return self._attributes

    @property
    def ids(self) -> Dict[str, NDArray]:
        """Dictionary of ids for the different shape types"""
        return self._ids

    def _update_ids(self):
        self._ids = {}
        for st in [Dot.name, Rectangle.name, Ellipse.name, Picture.name, PolygonShape.name]:
            self._ids[st] = np.flatnonzero(self._types == st)

    def todict(self) -> dict:
        """Dict representation of the shape array
        """

        return {"shape_array": [x.todict() for x in self.get_list()]}

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

        if np.any(self._types == PolygonShape.name):
            raise RuntimeError("tabular shape representation can not deal with "
                               "PolygonShapes")
        d = OrderedDict()
        d.update({"type": self._types,
                    "x": self._xy[:, 0].tolist(),
                "y": self._xy[:, 1].tolist(),
                "width": self._sizes[:, 0].tolist(),
                "height": self._sizes[:, 1].tolist(),
                "attributes": [str(x) for x in self._attributes.tolist()]
                })
        return d

    def add(self, shapes: Union[ShapeType, Tuple, Sequence, ShapeArray]):
        """add shape a the defined position

        Note
        ----
        see `add_somewhere` for adding shape at a random position
        """
        if isinstance(shapes, ShapeType):
            self._xy = np.append(self._xy, np.atleast_2d(shapes.xy), axis=0)
            self._sizes = np.append(
                self._sizes, np.atleast_2d(shapes.size), axis=0)
            self._attributes = np.append(self._attributes, shapes.attribute)
            self._types = np.append(self._types, shapes.name)
            self._polygons = np.append(self._polygons, shapes.polygon)
            self._convex_hull = None
            self._update_ids()

        elif isinstance(shapes, (list, tuple)):
            for x in shapes:
                self.add(x)

        elif isinstance(shapes, ShapeArray):
            self.add(shapes.get_list())

        else:
            raise TypeError(
                f"Can't add '{type(shapes)}', that's not a ShapeType or list of ShapeTypes."
            )

    def replace(self, index: int, shape: ShapeType):

        self._polygons[index] = shape.polygon
        self._xy[index, :] = shape.xy
        self._attributes[index] = shape.attribute
        self._sizes[index, :] = shape.size
        self._types[index] = shape.name
        self._convex_hull = None
        self._update_ids()

    def delete(self, index: IntOVector) -> None:
        self._polygons = np.delete(self._polygons, index)
        self._xy = np.delete(self._xy, index, axis=0)
        self._attributes = np.delete(self._attributes, index)
        self._types = np.delete(self._types, index)
        self._sizes = np.delete(self._sizes, index, axis=0)
        self._convex_hull = None
        self._update_ids()

    def clear(self):
        """ """
        self._polygons = np.empty(0, dtype=object)
        self._xy = np.empty((0, 2), dtype=np.float64)
        self._attributes = np.empty(0, dtype=object)
        self._types = np.empty(0, dtype=str)
        self._sizes = np.empty((0, 2), dtype=np.float64)
        self._convex_hull = None
        self._update_ids()

    def pop(self, index: Optional[int] = None) -> ShapeType:
        """Remove and return item at index"""
        if index is None:
            index = len(self._polygons) - 1
        rtn = self.get(index)
        self.delete(index)
        return rtn

    def sort_by_excentricity(self):
        """Sort order fo the objects in the array by excentricity, that is, by the
        distance the center of the convex hull (from close to far)
        """
        d = self.distances(Point2D(self.convex_hull.centroid))
        i = np.argsort(d)
        self._polygons = self._polygons[i]
        self._xy = self._xy[i, :]
        self._attributes = self._attributes[i]
        self._types = self._types[i]
        self._sizes = self._sizes[i]
        self._convex_hull = None
        self._update_ids()

    # def round_values(self, decimals: int = 0, int_type: type = np.int64,
    #                  rebuild_polygons=True) -> None: #FIXME rounding
    #     """rounds all values"""

    #     if decimals is None:
    #         return
    #     self._xy = round2(self._xy, decimals=decimals, int_type=int_type)
    #     self._dot_diameter = round2(self._dot_diameter, decimals=decimals,
    #                                 int_type=int_type)
    #     self._rect_sizes = round2(self._rect_sizes, decimals=decimals,
    #                               int_type=int_type)

    #     if rebuild_polygons:
    #         for i, shape in enumerate(self.get_list()):
    #             shape.delete_polygon()
    #             self._polygons[i] = shape.polygon

    def get(self, index: int) -> ShapeType:
        """Returns selected object"""
        name = self._types[index]
        if name == Dot.name:
            rtn = Dot(
                xy=self._xy[index, :],
                diameter=self._sizes[index, 0],
                attribute=self._attributes[index])
        elif name == Rectangle.name:
            return Rectangle(
                xy=self._xy[index, :],
                size=self._sizes[index, :],
                attribute=self._attributes[index])
        elif name == Picture.name:
            return Picture(
                xy=self._xy[index, :],
                size=self._sizes[index, :],
                path=self._attributes[index])  # type: ignore
        elif name == Ellipse.name:
            return Ellipse(
                xy=self._xy[index, :],
                size=self._sizes[index, :],
                attribute=self._attributes[index])
        elif name == PolygonShape.name:
            return PolygonShape(polygon=self._polygons[index],
                                attribute=self._attributes[index])
        else:
            raise TypeError(f"{name}: Unknown shape type id.")

        # rtn.set_polygon(self._polygons[index])
        return rtn

    def get_list(self, index: Optional[IntOVector] = None) -> List[ShapeType]:
        """Returns list with selected objects
        """
        if index is None:
            ids = range(self.n_objects)
        elif isinstance(index, int):
            ids = (index,)
        else:
            ids = index

        return [self.get(x) for x in ids]

    @property
    def n_objects(self) -> int:
        """number of shapes"""
        return len(self._attributes)

    @property
    def convex_hull(self) -> ConvexHull:
        """Convex hull poygon of the Shape Array"""
        if not isinstance(self._convex_hull, ConvexHull):
            self._convex_hull = ConvexHull(self._polygons)
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

        for i in range(self.n_objects):
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

    def distances(self, shape: Union[Point2D, ShapeType]) -> NDArray[np.float_]:
        """distances of a shape or Point2D to the elements in the shape array"""

        if isinstance(shape, (Point2D, CircularShapeType)):
            # circular target shape
            rtn = np.full(self.n_objects, np.nan)

            idx = self._ids[Dot.name]
            if len(idx) > 0:
                # circular -> dots in shape array
                rtn[idx] = _distance_circ_dot_array(obj=shape,
                                                    dots_xy=self._xy[idx, :],
                                                    dots_diameter=self._sizes[idx, 0])
            idx = self._ids[Ellipse.name]
            if len(idx) > 0:
                # circular -> ellipses in shape array
                rtn[idx] = _distance_circ_ellipse_array(obj=shape,
                                                        ellipses_xy=self._xy[idx, :],
                                                        ellipse_sizes=self._sizes[idx, :])
            # check if non-circular shapes are in shape_array
            idx = np.flatnonzero(np.isnan(rtn))
            if len(idx) > 0:
                if isinstance(shape, CircularShapeType):
                    rtn[idx] = shapely.distance(
                        shape.polygon, self._polygons[idx])
                else:
                    rtn[idx] = shapely.distance(
                        shape.xy_point, self._polygons[idx])
            return rtn

        else:
            # non-circular shape as target
            return shapely.distance(shape.polygon, self._polygons)

    def dwithin(self, shape: Union[Point2D, ShapeType],  distance: float = 0) -> NDArray[np.bool_]:
        """Returns True for all elements of the array that are within the
        specified distance.

        Note
        -----
        Using this function is more efficient than computing the distance and comparing the result.
        """
        if isinstance(shape, (Point2D, CircularShapeType)):
            rtn = np.full(self.n_objects, False)

            idx = self._ids[Dot.name]
            if len(idx) > 0:
                # circular -> dots in shape array
                dists = _distance_circ_dot_array(obj=shape,
                                                 dots_xy=self._xy[idx, :],
                                                 dots_diameter=self._sizes[idx, 0])
                rtn[idx] = dists < distance

            idx = self._ids[Ellipse.name]
            if len(idx) > 0:
                # circular -> ellipses in shape array
                dists = _distance_circ_ellipse_array(obj=shape,
                                                     ellipses_xy=self._xy[idx, :],
                                                     ellipse_sizes=self._sizes[idx, :])
                rtn[idx] = dists < distance

            # check if non-circular shapes are in shape_array
            idx = np.flatnonzero(np.isnan(rtn))
            if len(idx) > 0:
                if isinstance(shape, CircularShapeType):
                    rtn[idx] = shapely.dwithin(
                        shape.polygon, self._polygons[idx],  distance=distance)
                else:
                    rtn[idx] = shapely.dwithin(
                        shape.xy_point, self._polygons[idx],  distance=distance)

        else:
            # non-circular shape
            rtn = shapely.dwithin(
                shape.polygon, self._polygons, distance=distance)

        return rtn

    def contains_overlaps(self, min_distance: float = 0) -> bool:
        """Returns True for two or more elements overlap (i.e. taking
        into account the minimum distance).
        """
        for x in range(self.n_objects):
            if np.any(self.get_overlaps(x, min_distance)) > 0:
                return True

        return False

    def get_overlaps(self, index: int, min_distance: float = 0) -> NDArray[np.bool_]:
        """get overlaps with other objects. Ignores overlap with oneself."""
        overlaps = self.dwithin(self.get(index), distance=min_distance)
        overlaps[index] = False  # ignore overlap with oneself
        return overlaps

    def matrix_distances(self) -> NDArray:
        return self._relation_matrix(what=0)

    def matrix_dwithin(self, distance: float) -> NDArray:
        return self._relation_matrix(what=1, para=distance)

    def matrix_overlaps(self, min_distance: float = 0) -> NDArray:
        return self.matrix_dwithin(distance=min_distance)

    @staticmethod
    def from_dict(the_dict: Dict[str, Any]) -> ShapeArray:
        """read shape array from dict"""
        ##
        raise NotImplementedError()  # FIXME

    def _relation_matrix(self, what: int, para: float = 0) -> NDArray:
        """helper function returning the relation between polygons
        0 = distance
        1 = dwithin
        """
        arr = deepcopy(self)
        l = arr.n_objects
        rtn = np.full((l, l), np.nan)
        for x in reversed(range(l)):
            shape = arr.pop(x)
            if what == 0:
                y = arr.distances(shape)
            elif what == 1:
                y = arr.dwithin(shape=shape, distance=para)
            else:
                raise RuntimeError("unknown function")

            rtn[x, 0:x] = y

        # make symetric
        i_lower = np.triu_indices(l, 1)
        rtn[i_lower] = rtn.T[i_lower]
        return rtn

    def _random_free_position(self,
                              shape: ShapeType,
                              min_distance: float,
                              target_area: TargetArea,
                              ignore_overlaps: bool,
                              max_iterations: int) -> ShapeType:
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
                     max_iterations: int) -> int:  # FIXME manipulation objects
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

        target = self.get(index)

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


def _distance_circ_dot_array(obj: Union[Point2D, CircularShapeType],
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


def _distance_circ_ellipse_array(obj: Union[Point2D, CircularShapeType],
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
