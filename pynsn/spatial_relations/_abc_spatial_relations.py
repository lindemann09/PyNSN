__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from enum import Enum, auto
import numpy as np
from numpy.typing import NDArray

from .._lib.geometry import polar2cartesian
from pynsn._lib import geometry


class DisplTypes(Enum):
    """Displacement Types for spatial relation class """
    X = auto()
    Y = auto()
    RHO = auto()
    XY_SHORTEST = auto()
    SHORTEST = auto()


# all init-functions require 2D arrays
class ABCSpatialRelations(metaclass=ABCMeta):

    def __init__(self, a_xy: NDArray[np.floating],
                 b_xy: NDArray[np.floating],
                 a_relative_to_b: bool):
        """TODO """
        self._rho = None
        self._distances_rho = None
        self._distances_xy = None
        self._a_relative_to_b = a_relative_to_b
        if a_relative_to_b:
            self._xy_diff = np.atleast_2d(a_xy) - np.atleast_2d(b_xy)
        else:
            self._xy_diff = np.atleast_2d(b_xy) - np.atleast_2d(a_xy)

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
        """True if objects of array B are fully(!) inside the object of array A.
        if A_relative_to_B set to True, the function checks if arrays
        a objects of A are in B
        """

    @abstractmethod
    def displacement_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        """Distance of the required displacement of objects along the line
        between the object centers (rho) to have the minimum distance.

        Positive distance values indicate overlap and thus a required
        displacement in the direction rho to remove overlaps. Negative distance
        values indicate that the required displacement involves to move objects
        toward each other (not overlapping objects).
        """

    @property
    def is_above(self) -> NDArray:
        """Tests the relation of the object center. True if object center B is
        above object center A (or visa versa if A_relative_to_B =True)."""
        return self._xy_diff[:, 0] > 0

    @property
    def is_right(self) -> NDArray:
        """Tests the relation of the object center. True if object center B is
        right  of object center A (or visa versa if A_relative_to_B =True)."""
        return self._xy_diff[:, 1] > 0

    @property
    def rho(self) -> NDArray:
        """Polar coordinate rho of the line between the object centers. That is,
        position angles (in radians) of objects B relative to (viewed from) the
        objects A (or visa versa if A_relative_to_B =True).
        """
        if self._rho is None:
            self._rho = np.arctan2(self._xy_diff[:, 1],
                                   self._xy_diff[:, 0])
        return self._rho

    def displacements_polar(self, displ_type: DisplTypes,
                            minimum_gap: float = 0) -> NDArray:
        """Distance and angle of the required displacement of objects along the line
        between the object centers to have the minimum distance.

        Positive distance values indicate overlap and thus a required
        displacement in the direction rho to remove overlaps. Negative distance
        values indicate that the required displacement involves to move objects
        toward each other (not overlapping objects).

        use: required_displacement to get cartesian coordinates

        returns 2-d array (distance, angle)
        """
        # parent implements of DisplTypes.X and DisplTypes.Y
        if displ_type is DisplTypes.XY_SHORTEST:
            rtn = np.empty((len(self._xy_diff), 2, 2))
            rtn[:, :, 0] = self.displacements_polar(DisplTypes.X, minimum_gap)
            rtn[:, :, 1] = self.displacements_polar(DisplTypes.Y, minimum_gap)
            i = np.argmin(rtn[:, 0, :], axis=1)  # find min distances
            all_rows = np.arange(len(self._xy_diff))  # ":" does not work here
            return rtn[all_rows, :, i]
        elif displ_type is DisplTypes.SHORTEST:
            rtn = np.empty((len(self._xy_diff), 2, 3))
            rtn[:, :, 0] = self.displacements_polar(DisplTypes.X, minimum_gap)
            rtn[:, :, 1] = self.displacements_polar(DisplTypes.Y, minimum_gap)
            rtn[:, :, 2] = self.displacements_polar(
                DisplTypes.RHO, minimum_gap)
            i = np.argmin(rtn[:, 0, :], axis=1)  # find min distances
            all_rows = np.arange(len(self._xy_diff))  # ":" does not work here
            return rtn[all_rows, :, i]
        else:
            rtn = np.empty(self._xy_diff.shape)

            if displ_type is DisplTypes.X:
                rtn[:, 0] = -1*self.distances_xy[:, 0] + minimum_gap
                is_right = self._xy_diff[:, 0] > 0
                rtn[is_right, 1] = 0  # move to right
                rtn[~is_right, 1] = np.pi  # move to left
            elif displ_type is DisplTypes.Y:
                rtn[:, 0] = -1*self.distances_xy[:, 1] + minimum_gap
                is_above = self._xy_diff[:, 1] > 0
                rtn[is_above, 1] = geometry.NORTH  # move to top
                rtn[~is_above, 1] = geometry.SOUTH  # move to bottom
            elif displ_type is DisplTypes.RHO:
                # RHO
                rtn[:, 0] = self.displacement_distances_rho(
                    minimum_gap=minimum_gap)
                rtn[:, 1] = self.rho
            else:
                raise NotImplementedError()

            return rtn

    def displacements_cartesian(self, displ_type: DisplTypes,
                                minimum_gap: float = 0,
                                only_for_overlapping: bool = True) -> NDArray:
        """Required displacement of object B in cartesian coordinates
        to have no overlap with object A. (if A_relative_to_B = True,
        replacements of objects are return)

        Returns:
            replacements as array of vectors in cartesian coordinates
            if distance (replacements[:, 0]) is 0, no replacement required.

            Calculated movement direction (replacements[:, 0]) is always along the
            line between the two object center.
        """

        displ_polar = self.displacements_polar(displ_type=displ_type,
                                               minimum_gap=minimum_gap)
        if only_for_overlapping:
            # set distance =0, for all non override objects (dist<0)
            displ_polar[displ_polar[:, 0] < 0, 0] = 0

        return polar2cartesian(displ_polar)
