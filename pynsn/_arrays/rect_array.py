"""
Rectangle Array
"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from copy import deepcopy

from .._lib import misc
from .abc_object_array import ABCObjectArray
from .._lib.lib_typing import NumArray, OptArrayLike, IntOVector, Iterator, \
    Any, Union, Sequence, Optional, OptFloat, NDArray, ArrayLike
from .._shapes.rectangle import Rectangle
from .._shapes.dot import Dot
from .._lib.coordinate import Coordinate


class RectangleArray(ABCObjectArray):
    """Rectangle array class

    All object arrays are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 target_area: Union[Dot, Rectangle],
                 min_distance_between_objects: OptFloat = None,
                 min_distance_area_boarder: OptFloat = None,
                 xy: OptArrayLike = None,
                 sizes: OptArrayLike = None,
                 attributes: OptArrayLike = None) -> None:
        """Rectangular array is restricted to a certain area, it has a target area
        and a minimum gap.

        This properties allows shuffling free position and adapting
        properties.

        """
        super().__init__(xy=xy,
                         attributes=attributes,
                         target_area=target_area,
                         min_distance_between_objects=min_distance_between_objects,
                         min_distance_area_boarder=min_distance_area_boarder)
        self._sizes = np.array([])
        if sizes is not None:
            self._append_sizes(sizes)

        if self._xy.shape[0] != self._sizes.shape[0]:
            raise ValueError("Bad shaped data: " +
                             "xy has not the same length as sizes array")

    def _append_sizes(self, sizes: ArrayLike) -> int:
        """returns number of added rows"""
        sizes = misc.numpy_array_2d(sizes)
        if len(self._sizes) == 0:
            # ensure good shape of self.xy
            empty = np.array([]).reshape((0, 2))
            self._sizes = np.append(empty, sizes, axis=0)
        else:
            self._sizes = np.append(self._sizes, sizes, axis=0)
        return sizes.shape[0]

    def add(self, obj: Union[Rectangle, Sequence[Rectangle]]) -> None:
        """append one dot or list of dots"""
        if isinstance(obj, Rectangle):
            obj = [obj]
        else:
            obj = list(obj)

        for r in obj:
            assert isinstance(r, Rectangle)
            self._append_xy_attribute(xy=r.xy,
                                      attributes=r.attribute)
            self._append_sizes((r.width, r.height))
        self.properties.reset()

    @property
    def sizes(self) -> NDArray:
        """Numpy array of the object sizes

        The two dimensional array (shape=[2, n]) represents the width and height of the `n` objects in this array

        Returns:
            the width and height

        """
        return self._sizes

    @property
    def surface_areas(self) -> NDArray:
        # a = w*h
        return self._sizes[:, 0] * self._sizes[:, 1]

    @property
    def perimeter(self) -> NDArray:
        return 2 * (self._sizes[:, 0] + self._sizes[:, 1])

    def mod_round_values(self, decimals: int = 0,
                         int_type: type = np.int16) -> None:
        # inherited doc
        if decimals is None:
            return
        self._xy = misc.numpy_round2(self._xy, decimals=decimals,
                                     int_type=int_type)
        self._sizes = misc.numpy_round2(self._sizes, decimals=decimals,
                                        int_type=int_type)

    def to_dict(self) -> dict:
        # inherited doc
        d = super().to_dict()
        d.update({"sizes": self._sizes.tolist()})
        return d

    @staticmethod
    def from_dict(the_dict: dict) -> RectangleArray:
        """read rectangle array from dict"""

        rtn = RectangleArray(target_area=RectangleArray._target_area_from_dict(the_dict),
                             min_distance_between_objects=the_dict["min_distance_between_objects"],
                             min_distance_area_boarder=the_dict["min_distance_area_boarder"])
        rtn._append_xy_attribute(xy=the_dict["xy"],
                                 attributes=the_dict["attributes"])
        rtn._sizes = np.array(the_dict["sizes"])
        if len(rtn.sizes) != len(rtn.xy):
            raise RuntimeError("Badly shaped data: size data have not " +
                               "the same length as the coordinates")
        return rtn

    def clear(self) -> None:
        super().clear()
        self._sizes = np.array([])

    def delete(self, index: IntOVector) -> None:
        super().delete(index)
        self._sizes = np.delete(self._sizes, index, axis=0)

    def copy(self, indices: Union[int, Sequence[int], None] = None,
             deep_copy: bool = True) -> RectangleArray:
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._xy) == 0:
            return RectangleArray(
                target_area=deepcopy(self.target_area),
                min_distance_area_boarder=self.min_distance_area_boarder,
                min_distance_between_objects=self.min_distance_between_objects)
        if indices is None:
            indices = list(range(len(self._xy)))

        if deep_copy:
            return RectangleArray(
                target_area=deepcopy(self.target_area),
                min_distance_between_objects=self.min_distance_between_objects,
                min_distance_area_boarder=self.min_distance_area_boarder,
                xy=self._xy[indices, :].copy(),
                sizes=self._sizes[indices].copy(),
                attributes=self._attributes[indices].copy())
        else:
            return RectangleArray(
                target_area=deepcopy(self.target_area),
                min_distance_between_objects=self.min_distance_between_objects,
                min_distance_area_boarder=self.min_distance_area_boarder,
                xy=self._xy[indices, :],
                sizes=self._sizes[indices],
                attributes=self._attributes[indices])

    def _xy_distances(self, rect: Rectangle) -> NDArray:
        """return distances on both axes between rectangles and reference rec.
         negative number indicates overlap edges along that dimension.
        """
        if len(self._xy) == 0:
            return np.array([])
        else:
            pos_dist = np.abs(self._xy - rect.xy)
            max_not_overlap_dist = (self._sizes + rect.size) / 2
            dist = pos_dist - max_not_overlap_dist
            # FIXME intensive test distance function rect (also get_distance)
            return dist

    def get_distances(self, rect: Rectangle) -> NDArray:
        """Euclidean distances toward a single rectangle

        Negative numbers indicate an overlap of the objects

        Args:
            rect: object towards the distances will be calculated

        Returns: numpy array of distances
        """
        assert isinstance(rect, Rectangle)
        if len(self._xy) == 0:
            return np.array([])
        else:
            d_xy = self._xy_distances(rect)
            eucl_dist = np.hypot(d_xy[:, 0], d_xy[:, 1])
            for i, n_neg in enumerate(np.sum(d_xy < 0, axis=1)):
                if n_neg == 2:
                    # two dimensions overlap -> calc distance and make negative
                    eucl_dist[i] = -1 * eucl_dist[i]
                elif n_neg == 1:
                    # one dimension overlaps -> set to zero
                    if d_xy[i, 0] < 0:
                        x = 0
                    else:
                        x = d_xy[i, 0]
                    if d_xy[i, 1] < 0:
                        y = 0
                    else:
                        y = d_xy[i, 1]
                    eucl_dist[i] = np.hypot(x, y)

            return eucl_dist

    def iter_objects(self, indices: Optional[IntOVector] = None) -> Iterator[Rectangle]:
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
        if indices is None:
            data = zip(self._xy, self._sizes, self._attributes)
        else:
            data = zip(self._xy[indices, :], self._sizes[indices],
                       self._attributes[indices])

        for xy, s, att in data:
            rtn = Rectangle(xy=xy, size=s, attribute=att)
            yield rtn

    def find_objects(self, size: Optional[NumArray] = None,
                     attribute: Optional[Any] = None,
                     edge: Optional[Coordinate] = None) -> Sequence[int]:
        """returns indices of found objects

        2D-tuple
        """
        rtn = []
        for i in range(len(self._sizes)):
            if (size is not None and
                self._sizes[i, 0] != size[0] and self._sizes[i, 1] != size[1]) or \
                    (attribute is not None and self._attributes[i] != attribute):
                continue
            rtn.append(i)

        if edge is None:
            return rtn
        elif isinstance(edge, Coordinate):
            new_rtn = []
            for i, rect in zip(rtn, self.iter_objects(indices=rtn)):
                if edge in list(rect.iter_edges()):
                    new_rtn.append(i)
            return new_rtn
        else:
            raise TypeError("edge has to be a Coordinate and not {}.".format(
                type(edge)))

    def csv(self, variable_names: bool = True,
            hash_column: bool = False,
            attribute_column: bool = True) -> str:
        # inherited doc
        size_dict = {"width": self._sizes[:, 0], "height": self._sizes[:, 1]}
        if attribute_column:
            attr = self.attributes
        else:
            attr = None
        if hash_column:
            array_hash = self.hash
        else:
            array_hash = None
        return misc.make_csv(xy=self._xy,
                             size_data_dict=size_dict,
                             attributes=attr, array_hash=array_hash,
                             make_variable_names=variable_names)
