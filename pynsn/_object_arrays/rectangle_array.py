from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Any, Dict, Iterator, List, Optional, Sequence, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._shapes.coordinate import Coordinate
from .._lib.np_tools import make_vector_fixed_length, round2
from .._lib.typing import IntOVector
from .abc_object_aray import ABCObjectArray, hash_array, make_csv
from .._shapes.picture import Picture
from .._shapes.rectangle import Rectangle
from .._shapes import RectangleLike


class BaseRectangleArray:
    """Basic class comprising numpy lists of coordinates and sizes for rectangles"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None) -> None:
        if xy is not None and sizes is not None:
            self.xy = np.atleast_2d(xy)
            self.sizes = np.atleast_2d(sizes)
            if self.xy.shape != self.sizes.shape:
                raise ValueError("Badly shaped data: " +
                                 "xy has not the same shape as sizes array")
        else:
            self.xy = np.empty((0, 2))
            self.sizes = np.empty((0, 2))

    def n_objects(self) -> int:
        return len(self.xy)


class RectangleArray(BaseRectangleArray, ABCObjectArray):
    """Class of numpy representation of rectangles"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None,
                 attributes: Optional[ArrayLike] = None) -> None:
        super().__init__()
        self.attributes = np.array([])
        if xy is not None and sizes is not None:
            self.np_append(xy, sizes, attributes)

    def set_attributes(self, attributes: Optional[ArrayLike]) -> None:
        try:
            self.attributes = make_vector_fixed_length(
                attributes, length=self.xy.shape[0])
        except ValueError as err:
            raise ValueError("Length of attribute list does not match the " +
                             "size of the array.") from err

    def hash(self) -> str:
        return hash_array(xy=self.xy, perimeter=self.perimeter,
                          attributes=self.attributes)

    def np_append(self,
                  xy: ArrayLike,
                  sizes: ArrayLike,
                  attributes: Optional[ArrayLike] = None) -> None:
        """if attributes comprises only one element all new objects will get
        this attribute """
        xy = np.atleast_2d(xy)
        sizes = np.atleast_2d(sizes)
        if xy.shape != sizes.shape:
            raise ValueError("Badly shaped data: " +
                             "xy has not the same shape as sizes array")
        try:
            attributes = make_vector_fixed_length(attributes,
                                                  length=xy.shape[0])
        except ValueError as err:
            raise ValueError("Badly shaped data: attributes have not " +
                             "the correct length.") from err

        self.xy = np.append(self.xy, xy, axis=0)
        self.sizes = np.append(self.sizes, sizes, axis=0)
        self.attributes = np.append(self.attributes, attributes)

    def delete(self, index: IntOVector) -> None:
        self.xy = np.delete(self.xy, index, axis=0)
        self.sizes = np.delete(self.xy, index, axis=0)
        self.attributes = np.delete(self.attributes, index)

    def clear(self):
        self.xy = np.empty((0, 2))
        self.sizes = np.empty((0, 2))
        self.attributes = np.array([])

    def add(self, shapes: Union[RectangleLike, List[RectangleLike],
                                Sequence[RectangleLike]]) -> None:
        """append one rectangle/picture or list of rectangles/pictures"""
        if isinstance(shapes, Rectangle):
            shapes = [shapes]
        else:
            shapes = list(shapes)

        for obj in shapes:
            self.xy = np.append(self.xy, np.atleast_2d(obj.xy), axis=0)
            self.sizes = np.append(self.sizes, np.atleast_2d(obj.size), axis=0)
            self.attributes = np.append(self.attributes, obj.attribute)

    @ property
    def surface_areas(self) -> NDArray:
        # a = w*h
        return self.sizes[:, 0] * self.sizes[:, 1]

    @ property
    def perimeter(self) -> NDArray:
        return 2 * (self.sizes[:, 0] + self.sizes[:, 1])

    def round_values(self, decimals: int = 0,
                     int_type: type = np.int16) -> None:
        if decimals is None:
            return
        self.xy = round2(self.xy, decimals=decimals, int_type=int_type)
        self.sizes = round2(self.sizes, decimals=decimals, int_type=int_type)

    def iter(self, indices: Optional[IntOVector] = None) -> Iterator[RectangleLike]:
        """iterate over all or a part of the objects

        Parameters
        ----------
        indices

        Notes
        -----
        To iterate all object you might all use the class iterator __iter__:
        >>> for obj in my_array:
        >>>    print(obj)
        """

        if isinstance(indices, (int, np.integer)):
            fl_name = Picture.extract_filename(self.attributes[indices])
            if fl_name is None:
                yield Rectangle(xy=self.xy[indices, :],
                                size=self.sizes[indices],
                                attribute=self.attributes[indices])
            else:
                yield Picture(xy=self.xy[indices, :],
                              size=self.sizes[indices],
                              filename=fl_name)
        else:
            if indices is None:
                data = zip(self.xy,
                           self.sizes, self.attributes)
            else:
                data = zip(self.xy[indices, :],  # type: ignore
                           self.sizes[indices, :],  # type: ignore
                           self.attributes[indices])

            for xy, s, att in data:
                fl_name = Picture.extract_filename(att)
                if fl_name is None:
                    yield Rectangle(xy=xy, size=s, attribute=att)
                else:
                    yield Picture(xy=xy, size=s, filename=fl_name)

    def find(self,
             size: Union[Sequence[float],
                         NDArray[np.float_], None] = None,
             attribute: Optional[Any] = None,
             corner: Optional[Coordinate] = None) -> List[int]:
        """returns indices of found objects

        2D-tuple
        """
        rtn = []
        for i in range(len(self.sizes)):
            if (size is not None and
                self.sizes[i, 0] != size[0] and self.sizes[i, 1] != size[1]) or \
                    (attribute is not None and self.attributes[i] != attribute):
                continue
            rtn.append(i)

        if corner is None:
            return rtn
        elif isinstance(corner, Coordinate):
            new_rtn = []
            for i, rect in zip(rtn, self.iter(indices=rtn)):
                if corner in list(rect.iter_corners()):
                    new_rtn.append(i)
            return new_rtn
        else:
            raise TypeError("corner has to be a Coordinate and not {}.".format(
                type(corner)))

    def to_dict(self) -> dict:
        d = {"xy": self.xy.tolist()}
        d.update({"sizes": self.sizes.tolist()})
        if len(self.attributes) > 0:
            if np.all(self.attributes == self.attributes[0]):
                # all equal
                d.update({"attributes": self.attributes[0]})
            else:
                d.update({"attributes": self.attributes.tolist()})
        return d

    @classmethod
    def from_dict(cls, the_dict: Dict[str, Any]) -> RectangleArray:
        """read nsn stimulus from dict"""
        return cls(
            xy=the_dict["xy"],
            sizes=the_dict["sizes"],
            attributes=the_dict["attributes"])

    def copy(self,
             indices: Union[int, Sequence[int], NDArray[np.int_], None] = None,
             deep_copy: bool = True) -> RectangleArray:
        # inherited docs
        if indices is None or len(self.xy) == 0:
            xy = self.xy
            sizes = self.sizes
            attributes = self.attributes
        else:
            xy = self.xy[indices, :]
            sizes = self.sizes[indices]
            attributes = self.attributes[indices]

        if deep_copy:
            return RectangleArray(xy=xy.copy(),
                                  sizes=sizes.copy(),
                                  attributes=attributes.copy())
        else:
            return RectangleArray(xy=xy, sizes=sizes, attributes=attributes)

    def dataframe_dict(self,
                       hash_column: bool = False,
                       attribute_column: bool = True) -> dict:
        # inherited docs

        if hash_column:
            d = {"hash": [self.hash()] * len(self.xy)}
        else:
            d = {}
        d.update({"x": self.xy[:, 0].tolist(),
                  "y": self.xy[:, 1].tolist(),
                  "width": self.sizes[:, 0].tolist(),
                  "height": self.sizes[:, 1].tolist()})
        if attribute_column:
            d.update({"attributes": self.attributes.tolist()})
        return d

    def csv(self,
            variable_names: bool = True,
            hash_column: bool = False,
            attribute_column: bool = True) -> str:
        # inherited docs

        if attribute_column:
            attr = self.attributes
        else:
            attr = None
        if hash_column:
            array_hash = self.hash()
        else:
            array_hash = None

        return make_csv(xy=self.xy,
                        size_data_dict={"width": self.sizes[:, 0],
                                        "height": self.sizes[:, 1]},
                        attributes=attr,
                        array_hash=array_hash,
                        make_variable_names=variable_names)

    @property
    def center_of_mass(self) -> NDArray:
        """center of mass of all objects"""
        weighted_sum = np.sum(
            self.xy * np.atleast_2d(self.perimeter).T, axis=0)
        return weighted_sum / np.sum(self.perimeter)
