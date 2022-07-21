"""
Dot Array
"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np

from .abc_object_array import ABCObjectArray
from .. import misc
from ..lib_typing import OptInt, OptArrayLike, IntOVector, Iterator, \
    Any, Union, Sequence, Optional, OptFloat
from ..shapes.dot import Dot


# TODO: How to deal with rounding? Is saving to precises? Suggestion:
#  introduction precision parameter that is used by as_dict and get_csv and
#  hash


class DotArray(ABCObjectArray):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    def __init__(self,
                 target_area_radius: int,
                 min_dist_between: OptInt = None,
                 min_dist_area_boarder: OptInt = None,
                 xy: OptArrayLike = None,
                 diameter: OptArrayLike = None,
                 attributes: OptArrayLike = None) -> None:
        """Dot array is restricted to a certain area, it has a target area
        and a minimum gap.

        This properties allows shuffling free position and adapting
        properties.
        """
        super().__init__(xy=xy, attributes=attributes,
                         target_area_radius=target_area_radius,
                         min_dist_between=min_dist_between,
                         min_dist_area_boarder=min_dist_area_boarder)
        if diameter is None:
            self._diameter = np.array([])
        else:
            self._diameter = misc.numpy_vector(diameter)
        if self._xy.shape[0] != len(self._diameter):
            raise ValueError("Bad shaped data: " +
                             u"xy has not the same length as item_diameter")

    def add(self, dots: Union[Dot, Sequence[Dot]]) -> None:
        """append one dot or list of dots"""
        try:
            dots = list(dots)
        except TypeError:
            dots = [dots]
        for d in dots:
            assert isinstance(d, Dot)
            self._append_xy_attribute(xy=d.xy, attributes=d.attribute)
            self._diameter = np.append(self._diameter, d.diameter)
        self.properties.reset()

    @property
    def diameter(self) -> np.ndarray:
        return self._diameter

    @property
    def surface_areas(self) -> np.ndarray:
        # a = pi r**2 = pi d**2 / 4
        return np.pi * (self._diameter ** 2) / 4.0

    @property
    def perimeter(self) -> np.ndarray:
        return np.pi * self._diameter

    def mod_round_values(self, decimals: int = 0,
                         int_type=np.int32) -> None:
        """Round values of the array."""

        if decimals is None:
            return
        self._xy = misc.numpy_round2(self._xy, decimals=decimals,
                                     int_type=int_type)
        self._diameter = misc.numpy_round2(self._diameter, decimals=decimals,
                                           int_type=int_type)

    def as_dict(self) -> dict:
        """
        """
        d = super().as_dict()
        d.update({"diameter": self._diameter.tolist()})
        return d

    def read_from_dict(self, the_dict: dict) -> None:
        """read Dot collection from dict"""
        super().read_from_dict(the_dict)
        self._diameter = np.array(the_dict["diameter"])

    def clear(self) -> None:
        super().clear()
        self._diameter = np.array([])

    def delete(self, index: IntOVector) -> None:
        super().delete(index)
        self._diameter = np.delete(self._diameter, index)

    def copy(self, indices: OptArrayLike = None,
             deepcopy: bool = True) -> DotArray:
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._xy) == 0:
            return DotArray(target_area_radius=self.target_area_radius,
                            min_dist_between=self.min_dist_between,
                            min_dist_area_boarder=self.min_dist_area_boarder)

        if indices is None:
            indices = list(range(len(self._xy) ))

        if deepcopy:
            return DotArray(target_area_radius=self.target_area_radius,
                            min_dist_between=self.min_dist_between,
                            min_dist_area_boarder = self.min_dist_area_boarder,
                            xy=self._xy[indices, :].copy(),
                            diameter=self._diameter[indices].copy(),
                            attributes=self._attributes[indices].copy())
        else:
            return DotArray(target_area_radius=self.target_area_radius,
                            min_dist_between=self.min_dist_between,
                            min_dist_area_boarder = self.min_dist_area_boarder,
                            xy=self._xy[indices, :],
                            diameter=self._diameter[indices],
                            attributes=self._attributes[indices])

    def get_distances(self, dot: Dot) -> np.ndarray:
        """Distances toward a single dot
        negative numbers indicate overlap

        Returns
        -------
        distances : numpy array of distances
        """
        assert isinstance(dot, Dot)
        if len(self._xy) == 0:
            return np.array([])
        else:
            rtn = np.hypot(self._xy[:, 0] - dot.x, self._xy[:, 1] - dot.y) - \
                   ((self._diameter + dot.diameter) / 2.0)
            return rtn

    def iter_objects(self, indices: Optional[IntOVector] = None) -> Iterator[Dot]:
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
            yield Dot(xy=self._xy[indices, :],
                      diameter=self._diameter[indices],
                      attribute=self._attributes[indices])
        else:
            if indices is None:
                data = zip(self._xy, self._diameter, self._attributes)
            else:
                data = zip(self._xy[indices, :],  self._diameter[indices],
                                    self._attributes[indices])
            for xy, dia, att in data:
                yield Dot(xy=xy, diameter=dia, attribute=att)

    def find_objects(self, diameter: OptFloat = None,
                     attribute: Any = None) -> Sequence[int]:
        """returns indices of found objects
        """
        rtn = []
        for i in range(len(self._diameter)):
            if (diameter is not None and self._diameter[i] != diameter) or \
                    (attribute is not None and self._attributes[i] != attribute):
                continue
            rtn.append(i)
        return rtn

    def csv(self, variable_names: bool = True,
            hash_column: bool = False,
            attribute_column: bool = True) -> str:
        """Return the dot array as csv text

        Parameter
        ---------
        variable_names : bool, optional
            if True variable name will be printed in the first line
        """
        size_dict = {"diameter": self._diameter}
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