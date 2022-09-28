__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod

import numpy as np
from numpy.typing import NDArray

from .._lib.geometry import polar2cartesian

# all init-functions require 2D arrays


class ABCSpatialRelations(metaclass=ABCMeta):

    def __init__(self, a_xy: NDArray[np.floating],
                 b_xy: NDArray[np.floating],
                 reversed_relations: bool):
        self._radian = None
        self._distances = None
        self._reversed_relations = reversed_relations
        if reversed_relations:
            self._xy_diff = np.atleast_2d(b_xy) - np.atleast_2d(a_xy)
        else:
            self._xy_diff = np.atleast_2d(a_xy) - np.atleast_2d(b_xy)

    @property
    def angle(self) -> NDArray:
        """Angle (in radians) of objects A relative to the objects B
        (or visa versa if reversed_relations set to True)"""
        if self._radian is None:
            self._radian = np.arctan2(self._xy_diff[:, 1],
                                      self._xy_diff[:, 0])
        return self._radian

    @property
    @abstractmethod
    def distances(self) -> NDArray:
        """Euclidean distances between objects"""

    @abstractmethod
    def is_inside(self) -> NDArray:
        """True if objects of array B are fully(!) inside the object of array A.
        if reversed_relations set to True, the function checks if arrays
        a objects of A are in B
        """

    @abstractmethod
    def displacement_distances(self, minimum_distance: float = 0) -> NDArray:
        """Distance of the required displacement of objects along the line
        between the object centers to have the minimum distance.

        negative values indicate overlap and required displacement in the
        direction of the angle.

        use: required_displacement to get cartesian coordinates
        """

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        """Minimum required displacement of object B in cartesian coordinates
        to have no overlap with object A. (if reversed_relations = True,
        replacements of objects are return)

        Returns:
            replacements as array of vectors in cartesian coordinates
            if distance (replacements[:, 0]) is 0, no replacement required.

            Calculated movement direction (replacements[:, 0]) is always along the
            line between the two object center.
        """

        spatrel = np.array(
            (self.displacement_distances(minimum_distance=minimum_distance),
             self.angle)).T

        i = np.flatnonzero(spatrel[:, 0] < 0)  # displacement
        if self._reversed_relations:
            spatrel[i, 1] = np.pi + spatrel[i, 1]  # opposite direction

        # set all nan and override with overlapping relations
        rtn = np.zeros(spatrel.shape)
        rtn[i, :] = polar2cartesian(spatrel[i, :])
        return rtn
