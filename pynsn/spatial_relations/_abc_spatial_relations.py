__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from enum import Enum, auto

import numpy as np
from numpy.typing import NDArray
from .._lib import geometry


class SpreadTypes(Enum):
    """Displacement Types for spreading objects """
    X = auto()
    Y = auto()
    RHO = auto()
    CARDINAL = auto()  # shortest XY
    SHORTEST = auto()


class GatherTypes(Enum):
    """Displacement Types for gathering objects """
    RHO = auto()
    SHORTEST = auto()


# all init-functions require 2D arrays
class ABCSpatialRelations(metaclass=ABCMeta):

    def __init__(self,
                 a_xy: NDArray[np.floating],
                 b_xy: NDArray[np.floating],
                 a_relative_to_b: bool):
        """TODO """
        self._rho = None
        self._distances_rho = None
        self._distances_xy = None
        self._a_relative_to_b = a_relative_to_b
        if a_relative_to_b:
            self._xy_diff = np.atleast_2d(
                a_xy) - np.atleast_2d(b_xy)  # type: ignore
        else:
            self._xy_diff = np.atleast_2d(
                b_xy) - np.atleast_2d(a_xy)  # type: ignore

    @property
    @abstractmethod
    def distances_rho(self) -> NDArray:
        """Euclidean distances between objects along the line between the two
        object center."""

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
    def fits_inside(self, minimum_gap: float = 0) -> NDArray:
        """True if object B fits potentially inside object a (or visa versa)"""

    @abstractmethod
    def _spread_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        """Required displacement distance along the line between the object
        centers (rho) to have the minimum distance.

        Positive distance values indicate overlap and thus a required
        displacement in the direction rho to remove overlaps. Negative distance
        values indicate that the displacement would mean to move objects
        toward each other (i.e., two not overlapping objects).
        """

    @abstractmethod
    def _gather_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        """Required displacement distance along the line between the object
        centers (rho) to move object b inside object a.

        All distances are positive. Zero indicates no displacement required.
        NaN indicates that objects do not fit in each other.
        """

    @abstractmethod
    def _gather_polar_shortest(self, minimum_gap: float = 0) -> NDArray:
        """Polar coordinates of the shortest displacement that moves object B
        into object A (or visa versa)"""

    ### generic methods ###

    @property
    def is_above(self) -> NDArray:
        """Tests the relation of the object center. True if object center B is
        above object center A (or visa versa if a_relative_to_b =True)."""
        return self._xy_diff[:, 0] > 0

    @property
    def is_right(self) -> NDArray:
        """Tests the relation of the object center. True if object center B is
        right  of object center A (or visa versa if a_relative_to_b =True)."""
        return self._xy_diff[:, 1] > 0

    @property
    def rho(self) -> NDArray:
        """Polar coordinate rho of the line between the object centers. That is,
        position angles (in radians) of objects B relative to (viewed from) the
        objects A (or visa versa if a_relative_to_b =True).
        """
        if self._rho is None:
            self._rho = np.arctan2(self._xy_diff[:, 1],
                                   self._xy_diff[:, 0])
        return self._rho

    def spread(self,
               displ_type: SpreadTypes,
               minimum_gap: float = 0,
               polar: bool = False) -> NDArray:
        """The required displacement coordinates of objects to have the minimum
        distance.

        Positive distance values indicate overlap and thus a required
        displacement in the direction rho to remove overlaps. Negative distance
        values indicate that the required displacement involves to move objects
        toward each other (not overlapping objects).

        use: required_displacement to get cartesian coordinates

        returns 2-d array (distance, angle)
        """
        # parent implements of DisplTypes.X and DisplTypes.Y
        if displ_type is SpreadTypes.CARDINAL:
            rtn = np.empty((len(self._xy_diff), 2, 2))
            rtn[:, :, 0] = self.spread(SpreadTypes.X, minimum_gap, polar=True)
            rtn[:, :, 1] = self.spread(SpreadTypes.Y, minimum_gap, polar=True)
            i = np.argmin(rtn[:, 0, :], axis=1)  # find min distances
            all_rows = np.arange(len(self._xy_diff))  # ":" does not work here
            return rtn[all_rows, :, i]
        elif displ_type is SpreadTypes.SHORTEST:
            rtn = np.empty((len(self._xy_diff), 2, 3))
            rtn[:, :, 0] = self.spread(SpreadTypes.X, minimum_gap, polar=True)
            rtn[:, :, 1] = self.spread(SpreadTypes.Y, minimum_gap, polar=True)
            rtn[:, :, 2] = self.spread(
                SpreadTypes.RHO, minimum_gap, polar=True)
            i = np.argmin(rtn[:, 0, :], axis=1)  # find min distances
            all_rows = np.arange(len(self._xy_diff))  # ":" does not work here
            return rtn[all_rows, :, i]

        else:
            rtn = np.empty(self._xy_diff.shape)

            if displ_type is SpreadTypes.X:
                rtn[:, 0] = -1*self.distances_xy[:, 0] + minimum_gap
                is_right = self._xy_diff[:, 0] > 0
                rtn[is_right, 1] = 0  # move to right
                rtn[~is_right, 1] = np.pi  # move to left
            elif displ_type is SpreadTypes.Y:
                rtn[:, 0] = -1*self.distances_xy[:, 1] + minimum_gap
                is_above = self._xy_diff[:, 1] > 0
                rtn[is_above, 1] = geometry.NORTH  # move to top
                rtn[~is_above, 1] = geometry.SOUTH  # move to bottom
            elif displ_type is SpreadTypes.RHO:
                # RHO
                rtn[:, 0] = self._spread_distances_rho(
                    minimum_gap=minimum_gap)
                rtn[:, 1] = self.rho
            else:
                raise NotImplementedError()

            # remove non overlapping relations
            # set distance = 0, for all non override objects
            rtn[rtn[:, 0] < 0, 0] = 0
            if polar:
                return geometry.polar2cartesian(rtn)
            else:
                return rtn

    def gather(self,
               displ_type: GatherTypes,
               minimum_gap: float = 0,
               polar: bool = False) -> NDArray:
        """TODO"""

        if displ_type is GatherTypes.RHO:
            # RHO
            rtn = np.empty(self._xy_diff.shape)
            rtn[:, 0] = self._gather_distances_rho(
                minimum_gap=minimum_gap)
            rtn[:, 1] = self.rho

        elif displ_type is GatherTypes.SHORTEST:
            rtn = self._gather_polar_shortest(
                minimum_gap=minimum_gap)
        else:
            raise NotImplementedError()

        # remove non overlapping relations
        # set distance = 0, for all non override objects
        rtn[rtn[:, 0] < 0, 0] = 0
        if polar:
            return geometry.polar2cartesian(rtn)
        else:
            return rtn
