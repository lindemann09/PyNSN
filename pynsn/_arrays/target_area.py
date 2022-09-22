from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional

from .. import constants
from .._shapes import Dot, Rectangle, ShapeType


class TargetArea(object):

    def __init__(self,
                 target_area: ShapeType,
                 min_distance_between_objects: Optional[float] = None,
                 min_distance_area_boarder: Optional[float] = None) -> None:
        """TargetArea Parameter"""

        if not isinstance(target_area, (Dot, Rectangle)):
            raise ValueError("Area shape has to be a Dot or Rectangle")
        self.target_area = target_area
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
        if isinstance(self.target_area, Dot):
            area_size = self.target_area.diameter
        else:
            area_size = self.target_area.size
        return {"type": type(self).__name__,
                "target_area": type(self.target_area).__name__,
                "target_area_size": area_size,
                "min_distance_between_objects": self.min_distance_between_objects,
                "min_distance_area_boarder": self.min_distance_area_boarder}

    @staticmethod
    def _target_area_from_dict(dict_) -> ShapeType:
        # helper function
        if dict_["target_area"] == "Dot":
            return Dot(xy=(0, 0), diameter=dict_["target_area_size"])
        elif dict_["target_area"] == "Rectangle":
            return Rectangle(xy=(0, 0), size=dict_["target_area_size"])
        else:
            raise NotImplementedError()  # should never happen
