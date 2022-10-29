__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from typing import Union

import numpy as np
from numpy.typing import NDArray

from .._lib import geometry
from .._lib.exceptions import NoSolutionError
from .._object_arrays.dot_array import BaseDotArray
from .._object_arrays.rectangle_array import BaseRectangleArray


# all init-functions require 2D arrays
class ABCSpatialRelations(metaclass=ABCMeta):

    def __init__(self,
                 a_array: Union[BaseDotArray, BaseRectangleArray],
                 b_array: Union[BaseDotArray, BaseRectangleArray],
                 a_relative_to_b: bool):
        """TODO """
        self._rho = None
        self._distances_rho = None
        self._distances_xy = None
        self._a_relative_to_b = a_relative_to_b
        self._is_rectangle = (isinstance(a_array, BaseRectangleArray),
                              isinstance(b_array, BaseRectangleArray))
        if a_relative_to_b:
            self._xy_diff = np.atleast_2d(
                a_array.xy) - np.atleast_2d(b_array.xy)  # type: ignore
        else:
            self._xy_diff = np.atleast_2d(
                b_array.xy) - np.atleast_2d(a_array.xy)  # type: ignore

    @property
    @abstractmethod
    def distances(self) -> NDArray:
        """Shortest Euclidean distances between objects. Negative distances indicate
        overlap."""

    @property
    @abstractmethod
    def distances_rho(self) -> NDArray:
        """Euclidean distances between objects along the line between the two
        object center. Negative distances indicate overlap."""

    @property
    @abstractmethod
    def distances_xy(self) -> NDArray:
        """X an Y distances between objects."""

    @abstractmethod
    def is_inside(self) -> NDArray:
        """True if objects B are fully(!) inside the objects A.
        If a_relative_to_b set to True, the function checks if objects of A
        are in B.
        """

    @abstractmethod
    def fits_inside(self, minimum_gap: float) -> NDArray:
        """True if object B fits potentially inside object a (or visa versa)"""

    @abstractmethod
    def _spread_distances_rho(self, minimum_gap: float) -> NDArray:
        """Required displacement distance along the line between the object
        centers (rho) to have the minimum distance.

        Positive distance values indicate overlap and thus a required
        displacement in the direction rho to remove overlaps. Negative distance
        values indicate that the displacement would mean to move objects
        toward each other (i.e., two not overlapping objects).
        """

    @abstractmethod
    def _gather_displacements(self, minimum_gap: float, polar: bool) -> NDArray:
        """Cartesian coordinates of the shortest required displacement that
        moves object B into object A ( or visa versa)
        """

    ### generic methods ###

    @property
    def is_above(self) -> NDArray:
        """Tests the relation of the object center. True if object center B is
        above object center A ( or visa versa if a_relative_to_b =True)."""
        return self._xy_diff[:, 0] > 0

    @property
    def is_right(self) -> NDArray:
        """Tests the relation of the object center. True if object center B is
        right  of object center A ( or visa versa if a_relative_to_b =True)."""
        return self._xy_diff[:, 1] > 0

    @property
    def rho(self) -> NDArray:
        """Polar coordinate rho of the line between the object centers. That is ,
        position angles ( in radians) of objects B relative to (viewed from) the
        objects A ( or visa versa if a_relative_to_b =True).
        """
        if self._rho is None:
            self._rho = np.arctan2(self._xy_diff[:, 1],
                                   self._xy_diff[:, 0])
        return self._rho

    def gather(self, minimum_gap: float = 0, polar: bool = False) -> NDArray:
        """The required displacement coordinates of objects to moves object B
        into object A ( or visa versa)
        """
        if np.any(~self.fits_inside(minimum_gap)):
            n = np.sum(~self.fits_inside(minimum_gap))
            raise NoSolutionError("Not all objects can be gathered, "
                                  f"because some (n={n}) do not fit inside!")

        return self._gather_displacements(minimum_gap=minimum_gap, polar=polar)

    def spread(self,
               minimum_gap: float = 0,
               radial_displacements: bool = False,
               return_polar_coordinates: bool = False) -> NDArray:
        """The required displacement coordinates of objects to have the minimum
        distance.

        If objects move out of a
            * circular reference area: displacements will be radial, that is ,
              along the axes of the object center.
            * rectangular reference area: displacements will be the shortest
              displacements along the one of the cardinal axis(x or y), except
              `radial_displacements` is set to True.

        Positive distance values indicate overlap and thus a required
        displacement in the direction rho to remove overlaps. Negative distance
        values indicate that the required displacement involves to move objects
        toward each other(not overlapping objects).

        use: required_displacement to get cartesian coordinates

        returns 2-d array(distance, angle)
        """
        # FIXME A-relative-to-b not implemented
        # FIXME check logic of radial displacement

        if not radial_displacements and \
            ((not self._a_relative_to_b and self._is_rectangle[0]) or
             (self._a_relative_to_b and not self._is_rectangle[1])):
            # if objects moved out of rectangles area (that is, reference objects
            # are rectanglar) and if radial_displacement is not enforced
            # ----> cardinal displacements

            xy_move = np.empty((len(self._xy_diff), 2, 2))
            # x -> 0
            xy_move[:, 0, 0] = -1*self.distances_xy[:, 0] + minimum_gap
            is_right = self._xy_diff[:, 0] > 0
            xy_move[is_right, 1, 0] = 0  # move to right
            xy_move[~is_right, 1, 0] = np.pi  # move to left
            # y->1
            xy_move[:, 0, 1] = -1*self.distances_xy[:, 1] + minimum_gap
            is_above = self._xy_diff[:, 1] > 0
            xy_move[is_above, 1, 1] = geometry.NORTH  # move to top
            xy_move[~is_above, 1, 1] = geometry.SOUTH  # move to bottom

            # find shortest
            i = np.argmin(xy_move[:, 0, :], axis=1)  # find min distances
            all_rows = np.arange(len(self._xy_diff))  # ":" does not work here
            rtn_polar = xy_move[all_rows, :, i]

        else:
            # radial displacement
            rtn_polar = np.empty(self._xy_diff.shape)
            rtn_polar[:, 0] = self._spread_distances_rho(
                minimum_gap=minimum_gap)
            rtn_polar[:, 1] = self.rho

        # remove non overlapping relations
        # set distance = 0, for all non override objects
        rtn_polar[rtn_polar[:, 0] < 0, 0] = 0
        if not return_polar_coordinates:
            return geometry.polar2cartesian(rtn_polar)
        else:
            return rtn_polar
