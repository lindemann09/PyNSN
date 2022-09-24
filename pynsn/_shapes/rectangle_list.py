from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Any, Dict, Iterator, List, Optional, Sequence, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._lib.coordinate import Coordinate
from .._lib.np_tools import make_vector_fixed_length, round2
from .._lib.typing import IntOVector
from .abc_object_list import ABCObjectList
from .picture import Picture
from .rectangle import Rectangle

RectangleLike = Union[Rectangle, Picture]


class BaseRectangleList:
    """Basic class comprising numpy lists of coordinates and sizes for rectangles"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None) -> None:
        self.xy = np.empty((0, 2))
        self.sizes = np.empty((0, 2))

        if xy is not None and sizes is not None:
            self.xy = np.atleast_2d(xy)
            self.sizes = np.atleast_2d(sizes)
            self.check()

    def check(self):
        """raises value error if badly shaped data"""
        if self.xy.shape != self.sizes.shape:
            raise ValueError("Badly shaped data: " +
                             "xy has not the same shape as sizes array")


class RectangleList(BaseRectangleList, ABCObjectList):
    """Class of numpy representation of rectangles"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None,
                 attributes: Optional[ArrayLike] = None) -> None:
        super().__init__()
        self.attributes = np.array([])
        if xy is not None and sizes is not None:
            self.np_append(xy, sizes, attributes)

    def check(self) -> None:
        super().check()
        if len(self.attributes) != self.xy.shape[0]:
            raise ValueError("Badly shaped data: attributes have not " +
                             "the correct length.")

    def set_attributes(self, attributes: Optional[ArrayLike]) -> None:
        try:
            self.attributes = make_vector_fixed_length(
                attributes, length=self.xy.shape[0])
        except ValueError as err:
            raise ValueError("Length of attribute list does not match the " +
                             "size of the array.") from err

    def np_append(self,
                  xy: ArrayLike,
                  sizes: ArrayLike,
                  attributes: Optional[ArrayLike] = None) -> None:
        """if attributes comprises only one element all new objects will get
        this attribute """
        xy = np.atleast_2d(xy)
        try:
            attributes = make_vector_fixed_length(attributes,
                                                  length=xy.shape[0])
        except ValueError as err:
            raise ValueError("Badly shaped data: attributes have not " +
                             "the correct length.") from err

        self.xy = np.append(self.xy, xy, axis=0)
        self.sizes = np.append(self.sizes, np.atleast_2d(sizes), axis=0)
        self.attributes = np.append(self.attributes, attributes)
        self.check()

    def delete(self, index: IntOVector) -> None:
        self.xy = np.delete(self.xy, index, axis=0)
        self.sizes = np.delete(self.xy, index, axis=0)
        self.attributes = np.delete(self.attributes, index)

    def clear(self):
        self.xy = np.empty((0, 2))
        self.sizes = np.empty((0, 2))
        self.attributes = np.array([])

    def add(self, shapes: Union[Rectangle, Picture, List[RectangleLike],
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

    @staticmethod
    def from_dict(the_dict: Dict[str, Any]) -> RectangleList:
        """read nsn stimulus from dict"""
        return RectangleList(
            xy=the_dict["xy"],
            sizes=the_dict["sizes"],
            attributes=the_dict["attributes"])

    def copy(self,
             indices: Union[int, Sequence[int], NDArray[np.int_], None] = None,
             deep_copy: bool = True) -> RectangleList:
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
            return RectangleList(xy=xy.copy(),
                                 sizes=sizes.copy(),
                                 attributes=attributes.copy())
        else:
            return RectangleList(xy=xy, sizes=sizes, attributes=attributes)
