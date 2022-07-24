from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import json
from hashlib import md5

import numpy as np

from .._lib import geometry
from .._lib import misc
from .parameter import ArrayParameter
from .._lib.lib_typing import OptInt, OptArrayLike, ArrayLike, Union, \
    Sequence, Iterator, IntOVector, Optional
from .._shapes.dot import Dot
from .._shapes.point import Point
from .._shapes.rectangle import Rectangle
from .properties import ArrayProperties


class PointArray(ArrayParameter):
    """Class for attributes on two dimensional space"""

    def __init__(self,
                 target_area_radius: int,
                 min_dist_between: OptInt = None,
                 min_dist_area_boarder: OptInt = None,
                 xy: OptArrayLike = None,
                 attributes: OptArrayLike = None):
        """Numpy Position lists with attributes for optimized for numpy calculations

        Abstract class for implementation of dot and rect
        """
        super().__init__(target_area_radius=target_area_radius,
                         min_dist_between=min_dist_between,
                         min_dist_area_boarder=min_dist_area_boarder)

        self._xy = np.array([])
        self._attributes = np.array([])
        self._properties = ArrayProperties(self)

        if xy is not None:
            self._append_xy_attribute(xy=xy, attributes=attributes)

    def _append_xy_attribute(self, xy: ArrayLike,
                             attributes: OptArrayLike = None) -> int:
        """returns number of added rows"""
        xy = misc.numpy_array_2d(xy)
        if not isinstance(attributes, (tuple, list)):
            attributes = [attributes] * xy.shape[0]

        if len(attributes) != xy.shape[0]:
            raise RuntimeError("Badly shaped data: attributes have not " +
                               "the same length as the coordinates")

        self._attributes = np.append(self._attributes, attributes)
        if len(self._xy) == 0:
            empty = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
            self._xy = np.append(empty, xy, axis=0)
        else:
            self._xy = np.append(self._xy, xy, axis=0)
        self._properties.reset()
        return xy.shape[0]

    def add(self, points: Union[Point, Sequence[Point]]) -> None:
        """append one dot or list of dots"""
        try:
            points = list(points)
        except TypeError:
            points = [points]
        for p in points:
            assert isinstance(p, Point)
            self._append_xy_attribute(xy=p.xy, attributes=p.attribute)
        self.properties.reset()

    def __str__(self) -> str:
        prop_text = self._properties.as_text(extended_format=True)
        rtn = "- {}".format(type(self).__name__)
        rtn += "\n " + prop_text[1:]  # replace "-" with " "
        return rtn

    @property
    def xy(self) -> np.ndarray:
        return self._xy

    @property
    def xy_rounded_integer(self) -> np.ndarray:
        """rounded to integer"""
        return np.round(self._xy)

    @property
    def attributes(self) -> np.ndarray:
        return self._attributes

    @property
    def properties(self) -> ArrayProperties:
        return self._properties

    @property
    def surface_areas(self) -> np.ndarray:
        """per definition always zero"""
        return np.array([0] * len(self._xy))

    @property
    def perimeter(self) -> np.ndarray:
        """per definition always zero"""
        return np.array([0] * len(self._xy))

    def set_attributes(self, attributes: ArrayLike) -> None:
        """Set all attributes

        Parameter
        ---------
        attributes:  attribute (string) or list of attributes

        """

        if isinstance(attributes, (list, tuple)):
            if len(attributes) != self._properties.numerosity:
                raise ValueError("Length of attribute list does not adapt the " + \
                                 "size of the dot array.")
            self._attributes = np.array(attributes)
        else:
            self._attributes = np.array([attributes] * self._properties.numerosity)

    @property
    def hash(self) -> str:
        """md5_hash of positions and perimeter"""
        m = md5()
        m.update(
            self._xy.tobytes())  # to_byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        try:
            m.update(self.perimeter.tobytes())
        except AttributeError:
            pass
        m.update(self._attributes.tobytes())
        return m.hexdigest()

    def get_center_of_field_area(self) -> np.ndarray:
        """Center of all object positions
        """
        return geometry.center_of_positions(self.properties.convex_hull.xy)

    def clear(self) -> None:
        self._xy = np.array([])
        self._attributes = np.array([])
        self._properties.reset()

    def delete(self, index: ArrayLike) -> None:
        self._xy = np.delete(self._xy, index, axis=0)
        self._attributes = np.delete(self._attributes, index)
        self._properties.reset()

    def copy(self, indices: OptArrayLike = None,
             deepcopy: bool = True) -> PointArray:
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._xy) == 0:
            return PointArray(target_area_radius=self.target_area_radius,
                              min_dist_between=self.min_dist_between,
                              min_dist_area_boarder=self.min_dist_area_boarder)

        if indices is None:
            indices = list(range(len(self._xy)))

        if deepcopy:
            return PointArray(target_area_radius=self.target_area_radius,
                              min_dist_between=self.min_dist_between,
                              min_dist_area_boarder=self.min_dist_area_boarder,
                              xy=self._xy[indices, :].copy(),
                              attributes=self._attributes[indices].copy())
        else:
            return PointArray(target_area_radius=self.target_area_radius,
                              min_dist_between=self.min_dist_between,
                              min_dist_area_boarder=self.min_dist_area_boarder,
                              xy=self._xy[indices, :],
                              attributes=self._attributes[indices])

    def as_dict(self) -> dict:
        """
        """
        d = super().as_dict()
        d.update({"xy": self._xy.tolist()})
        if len(self._attributes) > 0 and misc.is_all_equal(self._attributes):
            d.update({"attributes": self._attributes[0]})
        else:
            d.update({"attributes": self._attributes.tolist()})
        return d

    @staticmethod
    def read_from_dict(the_dict: dict) -> PointArray:
        """read dot array from dict"""
        rtn = PointArray(target_area_radius=the_dict["target_area_radius"],
                         min_dist_between=the_dict["min_dist_between"],
                         min_dist_area_boarder=the_dict["min_dist_area_boarder"])
        rtn._append_xy_attribute(xy=the_dict["xy"],
                                 attributes=the_dict["attributes"])

        return rtn

    @staticmethod
    def load(json_file_name: str) -> PointArray:
        # override and extend read_from_dict not this function
        with open(json_file_name, 'r') as fl:
            return PointArray.read_from_dict(json.load(fl))


    def json(self, indent: int = None,
             include_hash: bool = False) -> str:
        """"""
        # override and extend as_dict not this function

        d = self.as_dict()
        if include_hash:
            d.update({"hash": self.hash})
        if not indent:
            indent = None
        return json.dumps(d, indent=indent)

    def save(self, json_file_name: str,
             indent: int = None,
             include_hash: bool = False) -> None:
        """"""
        with open(json_file_name, 'w') as fl:
            fl.write(self.json(indent=indent, include_hash=include_hash))

    def iter_objects(self, indices: Optional[IntOVector] = None) -> Iterator[Point]:
        """iterate over all or a part of the objects

        Parameters
        ----------
        indices: int or interable of integer

        Notes
        -----
        To iterate all object you might all use the class iterator __iter__:
        >>> for obj in my_array:
        >>>    print(obj)
        """

        if isinstance(indices, (int, np.integer)):
            yield Point(xy=self._xy[indices, :],
                      attribute=self._attributes[indices])
        else:
            if indices is None:
                data = zip(self._xy, self._attributes)
            else:
                data = zip(self._xy[indices, :],
                           self._attributes[indices])
            for xy, att in data:
                yield Point(xy=xy, attribute=att)

    def get_objects(self, indices: Sequence[int] = None) \
            -> Sequence[Union[Dot, Rectangle, Point]]:
        return list(self.iter_objects(indices=indices))

    def get_object(self, index: int) -> Union[None, Dot, Rectangle, Point]:
        if isinstance(index, int):
            return next(self.iter_objects(indices=index))
        else:
            raise ValueError("Index must be a integer not a {}. ".format(
                type(index).__name__) + "To handle multiple indices use 'get_objects'. ")

    def join(self, object_array) -> None:
        """add another object _arrays"""
        self.add(object_array.iter_objects())

