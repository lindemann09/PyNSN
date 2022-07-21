"""
Rectangle Array
"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np

from . import misc
from .base_classes import ABCObjectArray
from .lib_typing import OptArrayLike, IntOVector, ArrayLike, Iterator, \
    Any, Union, Sequence, Optional, NumPair
from .shapes import Rectangle
from .coordinate import Coordinate


class RectangleArray(ABCObjectArray):
    """
    """

    def __init__(self,
                 target_area_radius: int,
                 min_dist_between: bool = None,
                 min_dist_area_boarder: bool = None,
                 xy: OptArrayLike = None,
                 sizes: OptArrayLike = None,
                 attributes: OptArrayLike = None) -> None:
        """Rectangular array is restricted to a certain area, it has a target area
        and a minimum gap.

        This properties allows shuffling free position and adapting
        properties.

        """
        super().__init__(xy=xy, attributes=attributes,
                         target_area_radius=target_area_radius,
                         min_dist_between=min_dist_between,
                         min_dist_area_boarder=min_dist_area_boarder)
        self._sizes = np.array([])
        if sizes is not None:
            self._append_sizes(sizes)

        if self._xy.shape[0] != self._sizes.shape[0]:
            raise ValueError("Bad shaped data: " +
                             u"xy has not the same length as sizes array")

    def _append_sizes(self, sizes: ArrayLike) -> int:
        """returns number of added rows"""
        sizes = misc.numpy_array_2d(sizes)
        if len(self._sizes) == 0:
            empty = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
            self._sizes = np.append(empty, sizes, axis=0)
        else:
            self._sizes = np.append(self._sizes, sizes, axis=0)
        return sizes.shape[0]

    def add(self, rectangles: Union[Rectangle, Sequence[Rectangle]]) -> None:
        """append one dot or list of dots"""
        try:
            rectangles = list(rectangles)
        except TypeError:
            rectangles = [rectangles]

        for r in rectangles:
            assert isinstance(r, Rectangle)
            self._append_xy_attribute(xy=r.xy, attributes=r.attribute)
            self._append_sizes((r.width, r.height))
        self.properties.reset()

    @property
    def sizes(self) -> np.ndarray:
        return self._sizes

    @property
    def surface_areas(self) -> np.ndarray:
        # a = w*h
        return self._sizes[:, 0] * self._sizes[:, 1]

    @property
    def perimeter(self) -> np.ndarray:
        return 2 * (self._sizes[:, 0] + self._sizes[:, 1])

    def mod_round_values(self, decimals: int = 0,
                         int_type = np.int16) -> None:
        """Round values of the array."""

        if decimals is None:
            return
        self._xy = misc.numpy_round2(self._xy, decimals=decimals,
                                     int_type=int_type)
        self._sizes = misc.numpy_round2(self._sizes, decimals=decimals,
                                        int_type=int_type)

    def as_dict(self) -> dict:
        """
        """
        d = super().as_dict()
        d.update({"sizes": self._sizes.tolist()})
        return d

    def read_from_dict(self, the_dict: dict) -> None:
        """read rectangle array from dict"""
        super().read_from_dict(the_dict)
        self._sizes = np.array(the_dict["sizes"])

    def clear(self) -> None:
        super().clear()
        self._sizes = np.array([])

    def delete(self, index: IntOVector) -> None:
        super().delete(index)
        self._sizes = np.delete(self._sizes, index, axis=0)

    def copy(self, indices: OptArrayLike = None,
             deepcopy: bool = True) -> RectangleArray:
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._xy) == 0:
            return RectangleArray(
                target_area_radius=self.target_area_radius,
                min_dist_area_boarder=self.min_dist_area_boarder,
                min_dist_between=self.min_dist_between)
        if indices is None:
            indices = list(range(len(self._xy)))

        if deepcopy:
            return RectangleArray(
                target_area_radius=self.target_area_radius,
                min_dist_between=self.min_dist_between,
                min_dist_area_boarder=self.min_dist_area_boarder,
                xy=self._xy[indices, :].copy(),
                sizes=self._sizes[indices].copy(),
                attributes=self._attributes[indices].copy())
        else:
            return RectangleArray(
                target_area_radius=self.target_area_radius,
                min_dist_between=self.min_dist_between,
                min_dist_area_boarder=self.min_dist_area_boarder,
                xy=self._xy[indices, :],
                sizes=self._sizes[indices],
                attributes=self._attributes[indices])

    def _xy_distances(self, rect: Rectangle) -> np.ndarray:
        """return distances on both axes between rectangles and reference rec.
         negative number indicates overlap edges along that dimension.
        """
        if len(self._xy) == 0:
            return np.array([])
        else:
            pos_dist = np.abs(self._xy - rect.xy)
            max_not_overlap_dist = (self.sizes + rect.size) / 2
            dist = pos_dist - max_not_overlap_dist
            return dist # FIXME intensive test distance function rect (also get_distance)

    def get_distances(self, rect: Rectangle) -> np.ndarray:
        """Euclidean Distances toward a single Rectangle
        negative numbers indicate overlap

        Returns
        -------
        distances : numpy array of distances
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
            try:
                indices = list(indices)  # check if iterable
            except TypeError:
                indices = [indices]
            data = zip(self._xy[indices, :], self._sizes[indices],
                       self._attributes[indices])

        for xy, s, att in data:
            rtn = Rectangle(xy=xy, size=s, attribute=att)
            yield rtn

    def find_objects(self, size: Optional[NumPair] = None,
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
        """Return the rectangle array as csv text

        Parameter
        ---------
        variable_names : bool, optional
            if True variable name will be printed in the first line
        """
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
