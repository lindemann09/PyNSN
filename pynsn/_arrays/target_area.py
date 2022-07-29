from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from .._lib.lib_typing import OptInt
from .._lib import constants


class TargetArea(object):

    def __init__(self,
                 target_area_radius: int,
                 min_distance_between_objects: OptInt = None,
                 min_distance_area_boarder: OptInt = None) -> None:
        """TargetArea Parameter"""

        self.target_area_radius = target_area_radius
        if min_distance_between_objects is None:
            self.min_distance_between_objects = constants.DEFAULT_MIN_DISTANCE_BETWEEN_OBJECTS
        else:
            self.min_distance_between_objects = min_distance_between_objects
        if min_distance_area_boarder is None:
            self.min_distance_area_boarder = constants.DEFAULT_MIN_DISTANCE_AREA_BOARDER
        else:
            self.min_distance_area_boarder = min_distance_area_boarder

    def to_dict(self) -> dict:
        """Dict representation of the target area"""
        return {"type": type(self).__name__,
                "target_area_radius": self.target_area_radius,
                "min_distance_between_objects": self.min_distance_between_objects,
                "min_distance_area_boarder": self.min_distance_area_boarder}


# FIXME rectangular Target arrays
