"""

"""
from __future__ import annotations
from copy import deepcopy


__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .._shapes import (CircularShapeType, Dot, Ellipse, Picture, Point2D,
                       PolygonShape, Rectangle, ShapeType)
from .._shapes import ellipse_geometry as  ellipse_geo
from ..types import IntOVector
from .convex_hull import ConvexHull


class ShapeArray(object):
    """Numpy Optimizes Representation of Array of shapes"""

    SHAPE_LABELS = {
        Dot.ID: "Dot",
        Rectangle.ID: "Rectangle",
        Picture.ID: "Picture",
        Ellipse.ID: "Ellipse",
        PolygonShape.ID: "Polygon"}

    def __init__(self) -> None:
        self._xy = np.empty((0, 2), dtype=np.float64)
        self._sizes = np.empty((0, 2), dtype=np.float64)
        self._polygons = np.empty(0, dtype=object)
        self._attributes = np.empty(0, dtype=object)
        self._types = np.empty(0, dtype=int)
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
    def shape_types(self) -> List[str]:
        return [ShapeArray.SHAPE_LABELS[x] for x in self._types]

    @property
    def shape_types_ids(self) -> NDArray[np.int_]:
        return self._types

    @property
    def polygons(self) -> NDArray[shapely.Polygon]:
        return self._polygons

    @property
    def attributes(self) -> NDArray:
        return self._attributes

    @property
    def ids(self) -> Dict[int, NDArray]:
        """Dictionary of ids for the different shape types"""
        return self._ids

    def _update_ids(self):
        self._ids = {}
        for st in [Dot.ID, Rectangle.ID, Ellipse.ID, Picture.ID, PolygonShape.ID]:
            self._ids[st] = np.flatnonzero(self._types == st)

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
            self._types = np.append(self._types, shapes.ID)
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
        self._types[index] = shape.ID
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
        self._types = np.empty(0, dtype=int)
        self._sizes = np.empty((0, 2), dtype=np.float64)
        self._convex_hull = None
        self._update_ids()

    def pop(self, index: Optional[int]=None) -> ShapeType:
        """Remove and return item at index"""
        if index is None:
            index = len(self._polygons) -1
        rtn = self.get(index)
        self.delete(index)
        return rtn

    def sort_by_excentricity(self):
        ctr = self.convex_hull.centroid

        raise NotImplementedError()

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
        type_id = self._types[index]
        if type_id == Dot.ID:
            rtn = Dot(
                xy=self._xy[index, :],
                diameter=self._sizes[index, 0],
                attribute=self._attributes[index])
        elif type_id == Rectangle.ID:
            return Rectangle(
                xy=self._xy[index, :],
                size=self._sizes[index, :],
                attribute=self._attributes[index])
        elif type_id == Picture.ID:
            return Picture(
                xy=self._xy[index, :],
                size=self._sizes[index, :],
                path=self._attributes[index])  # type: ignore
        elif type_id == Ellipse.ID:
            return Ellipse(
                xy=self._xy[index, :],
                size=self._sizes[index, :],
                attribute=self._attributes[index])
        elif type_id == PolygonShape.ID:
            return PolygonShape(polygon=self._polygons[index],
                                attribute=self._attributes[index])
        else:
            raise TypeError(f"{type_id}: Unknown shape type id.")

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

    def to_dict(self) -> OrderedDict:
        """dict representation of the object array

        Notes. The can be used to create Pandas dataframe or Arrow Tables

        Examples
        --------
        >>> df_dict = stimulus.dataframe_dict()
        >>> df = pandas.DataFame(df_dict) # Pandas dataframe

        >>> table = pyarrow.Table.from_pydict(df_dict) # Arrow Table
        """

        d = OrderedDict()
        d.update({"x": self._xy[:, 0].tolist(),
                  "y": self._xy[:, 1].tolist(),
                  "width": self._sizes[:, 0].tolist(),
                  "height": self._sizes[:, 1].tolist(),
                  "attributes": self._attributes.tolist()})
        return d

    def csv(self, skip_columns: Optional[Sequence[str]] = None) -> str:
        """Comma-separated table representation the nsn stimulus

        Args:
            variable_names: if True first line include the variable names
            hash_column: if True hash will be included as first column
            skip_columns: list of column names that should not be exported
        Returns:
            CSV representation

        """

        d = self.to_dict()
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

            idx = self._ids[Dot.ID]
            if len(idx) > 0:
                # circular -> dots in shape array
                rtn[idx] = _distance_circ_dot_array(obj=shape,
                                                   dots_xy=self._xy[idx, :],
                                                   dots_diameter=self._sizes[idx, 0])
            idx = self._ids[Ellipse.ID]
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

    def dwithin(self, shape: Union[Point2D, ShapeType],  distance: float = 0) -> NDArray[np.int_]:
        """Returns True for all elements of the array that are within the
        specified distance.

        Note
        -----
        Using this function is more efficient than computing the distance and comparing the result.
        """
        if isinstance(shape, (Point2D, CircularShapeType)):
            rtn = np.full(self.n_objects, np.nan)

            idx = self._ids[Dot.ID]
            if len(idx) > 0:
                # circular -> dots in shape array
                dists = _distance_circ_dot_array(obj=shape,
                                                dots_xy=self._xy[idx, :],
                                                dots_diameter=self._sizes[idx, 0])
                rtn[idx] = dists < distance

            idx = self._ids[Ellipse.ID]
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


    def matrix_distances(self) -> NDArray:
        return _relation_matrix(self, what=0)

    def matrix_dwithin(self, distance:float) -> NDArray:
        return _relation_matrix(self, what=1, para=distance)


def from_dict(the_dict: Dict[str, Any]) -> ShapeArray:
    """read shape array from dict"""
    ##
    raise NotImplementedError() #FIXME


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
    theta = np.arctan2(d_xy[:, 0], d_xy[:, 1])
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


def _relation_matrix(shape_array: ShapeArray, what:int, para:float=0) -> NDArray:
    """helper function returning the relation between polygons
    0 = distance
    1 = dwithin
    """
    arr = deepcopy(shape_array)
    l = arr.n_objects
    rtn = np.full((l, l), np.nan)
    for x in reversed(range(l)):
        shape = arr.pop(x)
        if what==0:
            y = arr.distances(shape)
        elif what == 1:
            y = arr.dwithin(shape=shape, distance=para)
        else:
            raise RuntimeError("unkown function")

        rtn[x, 0:x] = y

    #make symetric
    i_lower = np.triu_indices(l, 1)
    rtn[i_lower] = rtn.T[i_lower]
    return rtn
