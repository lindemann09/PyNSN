"""

"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import json
from hashlib import md5

import numpy as np
from .._lib import geometry
from .._lib import misc
from .parameter import ArrayParameter
from .._lib.lib_typing import OptInt, OptArrayLike, Union, \
    Sequence, Iterator, IntOVector, Optional, NDArray, ArrayLike
from .._shapes.dot import Dot
from .._shapes.point import Point
from .._shapes.rectangle import Rectangle
from .properties import ArrayProperties

# FIXME PointArray not tested and not documented


class PointArray(ArrayParameter):
    """Point array class

    All object arrays are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 target_area_radius: int,
                 min_dist_between: OptInt = None,
                 min_dist_area_boarder: OptInt = None,
                 xy: OptArrayLike = None,
                 attributes: OptArrayLike = None) -> None:
        # also as parent  class for implementation of dot and rect arrays

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
            attributes = np.array([attributes] * xy.shape[0])

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

    def add(self, obj: Union[Point, Sequence[Point]]) -> None:
        """append one dot or list of dots"""
        if isinstance(obj, Point):
            obj = [obj]
        else:
            obj = list(obj)

        for p in obj:
            assert isinstance(p, Point)
            self._append_xy_attribute(xy=p.xy, attributes=p.attribute)
        self.properties.reset()

    def __str__(self) -> str:
        prop_text = self._properties.as_text(extended_format=True)
        rtn = "- {}".format(type(self).__name__)
        rtn += "\n " + prop_text[1:]  # replace "-" with " "
        return rtn

    @property
    def xy(self) -> NDArray:
        """Numpy array of the object locations

        The two dimensional array (shape=[2, `n`]) represents the locations of the center of
        the `n` objects in this array
        """
        return self._xy

    @property
    def attributes(self) -> NDArray:
        """Numpy vector of the object attributes
        """
        return self._attributes

    @property
    def properties(self) -> ArrayProperties:
        """Properties of the object array.

        ``ArrayProperties`` represents and handles (fitting, scaling) visual
        properties of the object like

        * numerosity
        * average_dot_diameter/average_rectangle_size
        * total_surface_area
        * average_surface_area
        * total_perimeter
        * average_perimeter
        * field_area
        * field_area_positions
        * sparsity
        * log_spacing
        * log_size
        * coverage
        """
        return self._properties

    @property
    def surface_areas(self) -> NDArray:
        """Size of all points is per definition always zero"""
        return np.array([0] * len(self._xy))

    @property
    def perimeter(self) -> NDArray:
        """Perimeter of all points is per definition always zero"""
        return np.array([0] * len(self._xy))

    def set_attributes(self, attributes: ArrayLike) -> None:
        """Set all attributes

        Args:
            attributes: attribute (string) or list of attributes
        """
        if isinstance(attributes, (list, tuple)):
            if len(attributes) != self._properties.numerosity:
                raise ValueError("Length of attribute list does not adapt the " +
                                 "size of the dot array.")
            self._attributes = np.array(attributes)
        else:
            self._attributes = np.array([attributes] * self._properties.numerosity)

    @property
    def hash(self) -> str:
        """Hash (MD5 hash) of the array

        The hash can be used as an unique identifier of the object array.

        Notes:
            Hashing is based on the byte representations of the positions, perimeter
            and attributes.
        """
        m = md5()
        m.update(
            self._xy.tobytes())  # to_byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        try:
            m.update(self.perimeter.tobytes())
        except AttributeError:
            pass
        m.update(self.attributes.tobytes())
        return m.hexdigest()

    def get_center_of_field_area(self) -> NDArray:
        """Center of all objects

        Returns:
            Coordinate of center position
        """
        return geometry.center_of_positions(self.properties.convex_hull.xy)

    def clear(self) -> None:
        self._xy = np.array([])
        self._attributes = np.array([])
        self._properties.reset()

    def delete(self, index: IntOVector) -> None:
        """Delete object

        Args:
            index:

        Returns:

        """
        self._xy = np.delete(self._xy, index, axis=0)
        self._attributes = np.delete(self._attributes, index)
        self._properties.reset()

    def copy(self, indices: Union[int, Sequence[int], None] = None,
             deep_copy: bool = True) -> PointArray:
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._xy) == 0:
            return PointArray(target_area_radius=self.target_area_radius,
                              min_dist_between=self.min_dist_between,
                              min_dist_area_boarder=self.min_dist_area_boarder)

        if indices is None:
            indices = list(range(len(self._xy)))

        if deep_copy:
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

    def to_dict(self) -> dict:
        """Convert array tp dictionary
        """
        d = super().to_dict()
        d.update({"xy": self._xy.tolist()})
        if len(self._attributes) > 0 and misc.is_all_equal(self._attributes):
            d.update({"attributes": self._attributes[0]})
        else:
            d.update({"attributes": self._attributes.tolist()})
        return d

    @staticmethod
    def from_dict(the_dict: dict) -> PointArray:
        """read dot array from dict"""
        rtn = PointArray(target_area_radius=the_dict["target_area_radius"],
                         min_dist_between=the_dict["min_dist_between"],
                         min_dist_area_boarder=the_dict["min_dist_area_boarder"])
        rtn._append_xy_attribute(xy=the_dict["xy"],
                                 attributes=the_dict["attributes"])

        return rtn

    def json(self, indent: Optional[int] = None,
             include_hash: bool = False) -> str:
        """"""
        # override and extend to_dict not this function

        d = self.to_dict()
        if include_hash:
            d.update({"hash": self.hash})
        if not indent:
            indent = None
        return json.dumps(d, indent=indent)

    def save(self, json_file_name: str,
             indent: Optional[int] = None,
             include_hash: bool = False) -> None:
        """"""
        with open(json_file_name, 'w') as fl:
            fl.write(self.json(indent=indent, include_hash=include_hash))

    def iter_objects(self, indices: Optional[IntOVector] = None) -> Iterator[Point]:
        """iterate over all or a part of the objects

        Parameters
        ----------
        indices: int or iterable of integer

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

    def get_objects(self, indices: Optional[Sequence[int]] = None) \
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
