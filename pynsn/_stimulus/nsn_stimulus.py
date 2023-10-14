"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from copy import deepcopy
import warnings
from typing import Union

from .shapes import Dot, Rectangle
from .properties import ArrayProperties
from .shape_array import ShapeArray

from .._lib.misc import key_value_format


class NSNStimulus(object):
    """Non-Symbolic Number Stimulus

    NSN-Stimulus are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(
        self,
        target_area: Union[Dot, Rectangle],
        min_dist_between_shapes: float = 2,
        min_dist_area_edge: float = 2,
    ) -> None:
        super().__init__()
        assert isinstance(target_area, (Dot, Rectangle))

        if tuple(target_area.xy) != (0, 0):
            warnings.warn("TargetArea does not use shape position. "
                          "Shape Position will be set to (0, 0).",
                          UserWarning)

        self.target_area = target_area.variant(xy=(0, 0))
        self.min_dist_between_objects = min_dist_between_shapes
        self.min_dist_area_edge = min_dist_area_edge
        self._shapes = ShapeArray()

    @property
    def shapes(self) -> ShapeArray:
        return self._shapes

    @property
    def properties(self) -> ArrayProperties:
        return self._shapes.properties

    def deepcopy(self) -> NSNStimulus:
        """A (deep) copy of the nsn stimulus.

        Returns:
            a copy of the array
        """
        rtn = NSNStimulus(
            target_area=deepcopy(self.target_area),
            min_dist_between_shapes=self.min_dist_area_edge,
            min_dist_area_edge=self.min_dist_area_edge
        )
        rtn.shapes.add(self._shapes)
        return rtn

    def properties_txt(self, with_hash: bool = False, extended_format: bool = False) -> str:
        prop = self.properties
        if extended_format:
            rtn = None
            for k, v in prop.to_dict().items():
                if rtn is None:
                    if with_hash:
                        rtn = f"-{k}  {v}\n "
                    else:
                        rtn = "-"
                else:
                    rtn += key_value_format(k, v) + "\n "

            if rtn is None:
                rtn = ""
        else:
            if with_hash:
                rtn = "ID: {} ".format(self._shapes.hash())
            else:
                rtn = ""
            rtn += (
                f"N: {prop.numerosity}, "
                + f"TSA: {int(prop.total_surface_area)}, "
                + f"ISA: {int(prop.average_surface_area)}, "
                + f"FA: {int(prop.field_area)}, "
                + f"SPAR: {prop.sparsity:.2f}, "
                + f"logSIZE: {prop.log_size:.2f}, "
                + f"logSPACE: {prop.log_spacing:.2f}, "
                + f"COV: {prop.coverage:.2f}"
            )

        return rtn.rstrip()
