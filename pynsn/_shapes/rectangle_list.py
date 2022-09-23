__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional, Sequence, Union, Iterator

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._lib.typing import IntOVector
from .._lib.np_tools import make_vector_fixed_length, round2

from .rectangle import Rectangle
from .picture import Picture


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


class RectangleList(BaseRectangleList):
    """Class of numpy representation of rectangles"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None,
                 attributes: Optional[ArrayLike] = None) -> None:
        super().__init__()
        self.attributes = np.array([])
        if xy is not None and sizes is not None:
            self.append(xy, sizes, attributes)

    def check(self) -> None:
        """raises value error if badly shaped data"""
        super().check()
        if len(self.attributes) != self.xy.shape[0]:
            raise ValueError("Badly shaped data: attributes have not " +
                             "the correct length.")

    def set_attributes(self, attributes: Optional[ArrayLike]) -> None:
        """Set all attributes

        Args:
            attributes: single attribute or list of attributes
        """
        try:
            self.attributes = make_vector_fixed_length(
                attributes, length=self.xy.shape[0])
        except ValueError as err:
            raise ValueError("Length of attribute list does not match the " +
                             "size of the array.") from err

    def append(self,
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

    def add_shapes(self, shapes: Union[Rectangle, Picture,
                                       Sequence[Union[Rectangle, Picture]]]) -> None:
        """append one dot or list of dots"""
        if isinstance(shapes, Rectangle):
            shapes = [shapes]
        else:
            shapes = list(shapes)

        for obj in shapes:
            self.xy = np.append(self.xy, np.atleast_2d(obj.xy), axis=0)
            self.sizes = np.append(self.sizes, np.atleast_2d(obj.size), axis=0)
            self.attributes = np.append(self.attributes, obj.attribute)

    @property
    def surface_areas(self) -> NDArray:
        # a = w*h
        return self.sizes[:, 0] * self.sizes[:, 1]

    @property
    def perimeter(self) -> NDArray:
        return 2 * (self.sizes[:, 0] + self.sizes[:, 1])

    def round_values(self, decimals: int = 0,
                     int_type: type = np.int16) -> None:
        """Round all values

        Args:
            decimals: number of decimal places
            int_type: numpy int type (default=numpy.int16)
        """
        if decimals is None:
            return
        self.xy = round2(self.xy, decimals=decimals, int_type=int_type)
        self.sizes = round2(self.sizes, decimals=decimals, int_type=int_type)

    def iter_objects(self, indices: Optional[IntOVector] = None) \
            -> Iterator[Union[Rectangle, Picture]]:
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
