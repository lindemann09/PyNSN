from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from ..lib_typing import OptInt
from ... import constants


class ArrayParameter(object):

    def __init__(self,
                 target_area_radius: int,
                 min_dist_between: OptInt = None,
                 min_dist_area_boarder: OptInt = None) -> None:
        """Numpy Position lists with attributes for optimized for numpy calculations

        Abstract class for implementation of dot and rect
        """
        self.target_area_radius = target_area_radius
        if min_dist_between is None:
            self.min_dist_between = constants.DEFAULT_MIN_DIST_BETWEEN
        else:
            self.min_dist_between = min_dist_between
        if min_dist_area_boarder is None:
            self.min_dist_area_boarder = constants.DEFAULT_MIN_DIST_AREA_BOARDER
        else:
            self.min_dist_area_boarder = min_dist_area_boarder

    def as_dict(self) -> dict:
        return {"type": type(self).__name__,
                "target_area_radius": self.target_area_radius,
                "min_dist_between": self.min_dist_between,
                "min_dist_area_boarder": self.min_dist_area_boarder}

