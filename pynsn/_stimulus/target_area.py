from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional

from .. import constants
from .._shapes import Dot, Rectangle, ShapeType


class TargetArea(object):

    def __init__(self,
                 target_area_shape: ShapeType,
                 min_dist_between_objects: Optional[float] = None,
                 min_dist_area_edge: Optional[float] = None) -> None:
        """TargetArea Parameter"""

        if not isinstance(target_area_shape, (Dot, Rectangle)):
            raise ValueError("Area shape has to be a Dot or Rectangle")
        self.target_area_shape = target_area_shape
        if min_dist_between_objects is None:
            self.min_dist_between_objects = constants.DEFAULT_min_dist_between_objects
        else:
            self.min_dist_between_objects = min_dist_between_objects
        if min_dist_area_edge is None:
            self.min_dist_area_edge = constants.DEFAULT_min_dist_area_edge
        else:
            self.min_dist_area_edge = min_dist_area_edge

    def to_dict(self) -> dict:
        """Dict representation of the target area"""
        if isinstance(self.target_area_shape, Dot):
            area_size = self.target_area_shape.diameter
        else:
            area_size = self.target_area_shape.size
        return {"target_area_shape": type(self.target_area_shape).__name__,
                "target_area_size": area_size,
                "min_dist_between_objects": self.min_dist_between_objects,
                "min_dist_area_edge": self.min_dist_area_edge}

    @classmethod
    def from_dict(cls, the_dict: dict) -> TargetArea:
        if the_dict["target_area_shape"] == "Dot":
            shape = Dot(xy=(0, 0), diameter=the_dict["target_area_size"])
        elif the_dict["target_area_shape"] == "Rectangle":
            shape = Rectangle(xy=(0, 0), size=the_dict["target_area_size"])
        else:
            raise NotImplementedError()  # should never happen
        return cls(
            target_area_shape=shape,
            min_dist_area_edge=the_dict["min_dist_area_edge"],
            min_dist_between_objects=the_dict["min_dist_between_objects"])
