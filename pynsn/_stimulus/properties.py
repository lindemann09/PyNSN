"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import enum
from collections import OrderedDict
from typing import Any, Union

import numpy as np
from numpy.typing import NDArray

from .. import _stimulus
from .convex_hull import ConvexHull


class VisualPropertyFlags(enum.IntFlag):
    AV_DOT_DIAMETER = enum.auto()
    AV_SURFACE_AREA = enum.auto()
    AV_PERIMETER = enum.auto()
    AV_RECT_SIZE = enum.auto()

    TOTAL_SURFACE_AREA = enum.auto()
    TOTAL_PERIMETER = enum.auto()
    SPARSITY = enum.auto()
    FIELD_AREA = enum.auto()
    COVERAGE = enum.auto()

    LOG_SPACING = enum.auto()
    LOG_SIZE = enum.auto()

    NUMEROSITY = enum.auto()

    def is_dependent_from(self, other_property: Any) -> bool:
        """returns true if both properties are not independent"""
        return (self.is_size_property() and other_property.is_size_property()) or \
               (self.is_space_property() and other_property.is_space_property())

    def is_size_property(self) -> bool:
        return self in (VisualPropertyFlags.LOG_SIZE,
                        VisualPropertyFlags.TOTAL_SURFACE_AREA,
                        VisualPropertyFlags.AV_DOT_DIAMETER,
                        VisualPropertyFlags.AV_SURFACE_AREA,
                        VisualPropertyFlags.AV_PERIMETER,
                        VisualPropertyFlags.TOTAL_PERIMETER)

    def is_space_property(self) -> bool:
        return self in (VisualPropertyFlags.LOG_SPACING,
                        VisualPropertyFlags.SPARSITY,
                        VisualPropertyFlags.FIELD_AREA)

    def label(self) -> str:
        labels = {
            VisualPropertyFlags.NUMEROSITY: "Numerosity",
            VisualPropertyFlags.LOG_SIZE: "Log size",
            VisualPropertyFlags.TOTAL_SURFACE_AREA: "Total surface area",
            VisualPropertyFlags.AV_DOT_DIAMETER: "Av. dot diameter",
            VisualPropertyFlags.AV_SURFACE_AREA: "Av. surface area",
            VisualPropertyFlags.AV_PERIMETER: "Av. perimeter",
            VisualPropertyFlags.TOTAL_PERIMETER: "Total perimeter",
            VisualPropertyFlags.AV_RECT_SIZE: "Av. rectangle Size",
            VisualPropertyFlags.LOG_SPACING: "Log spacing",
            VisualPropertyFlags.SPARSITY: "Sparsity",
            VisualPropertyFlags.FIELD_AREA: "Field area",
            VisualPropertyFlags.COVERAGE: "Coverage"}
        return labels[self]


class ArrayProperties(object):
    """Non-Symbolic Number Stimulus"""

    def __init__(self, shape_array: _stimulus.ShapeArray) -> None:
        self._shapes = shape_array
        self._ch = None

    @property
    def numerosity(self) -> int:
        """number of shapes"""
        return self._shapes.n_objects

    def reset_convex_hull(self) -> None:
        """reset to enforce recalculation of convex hull"""
        self._ch = None  # convex_hull

    @property
    def convex_hull(self) -> ConvexHull:
        if not isinstance(self._ch, ConvexHull):
            self._ch = ConvexHull(self._shapes.polygons)
        return self._ch

    @property
    def average_dot_diameter(self) -> float:
        rtn = np.nanmean(self._shapes.dot_diameter)
        if np.isnan(rtn):
            return 0.0
        return float(rtn)

    @property
    def average_rectangle_size(self) -> NDArray:
        rtn = np.nanmean(self._shapes.rect_sizes, axis=0)
        if np.isnan(rtn[0]):
            return np.array([0.0, 0.0])
        return rtn

    @property
    def total_surface_area(self) -> float:
        return np.nansum(self._shapes.areas)

    @property
    def average_surface_area(self) -> float:
        if self._shapes.n_objects == 0:
            return 0.0
        return np.nanmean(self._shapes.areas)

    @property
    def total_perimeter(self) -> float:
        return np.nansum(self._shapes.perimeter)

    @property
    def average_perimeter(self) -> float:
        if self._shapes.n_objects == 0:
            return 0.0
        return np.nanmean(self._shapes.perimeter)

    @property
    def coverage(self) -> Union[np.floating, float]:
        """percent coverage in the field area. It takes thus the object size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity
        """
        try:
            return self.total_surface_area / self.field_area
        except ZeroDivisionError:
            return np.nan

    @property
    def log_size(self) -> float:
        try:
            return np.log2(self.total_surface_area) + np.log2(self.average_surface_area)
        except ValueError:
            return np.nan

    @property
    def log_spacing(self) -> float:
        try:
            return np.log2(self.field_area) + np.log2(self.sparsity)
        except ValueError:
            return np.nan

    @property
    def sparsity(self) -> float:
        try:
            return self.field_area / self.numerosity
        except ZeroDivisionError:
            return np.nan

    @property
    def field_area(self) -> float:
        return self.convex_hull.area

    # def get(self, property_flag: VisualPropertyFlags) -> Any:
    #     """returns a visual property"""

    #     assert isinstance(property_flag, VisualPropertyFlags)

    #     # Adapt
    #     if property_flag == VisualPropertyFlags.AV_DOT_DIAMETER:
    #         return self.average_dot_diameter

    #     elif property_flag == VisualPropertyFlags.AV_RECT_SIZE:
    #         return self.average_rectangle_size

    #     elif property_flag == VisualPropertyFlags.AV_PERIMETER:
    #         return self.average_perimeter

    #     elif property_flag == VisualPropertyFlags.TOTAL_PERIMETER:
    #         return self.total_perimeter

    #     elif property_flag == VisualPropertyFlags.AV_SURFACE_AREA:
    #         return self.average_surface_area

    #     elif property_flag == VisualPropertyFlags.TOTAL_SURFACE_AREA:
    #         return self.total_surface_area

    #     elif property_flag == VisualPropertyFlags.LOG_SIZE:
    #         return self.log_size

    #     elif property_flag == VisualPropertyFlags.LOG_SPACING:
    #         return self.log_spacing

    #     elif property_flag == VisualPropertyFlags.SPARSITY:
    #         return self.sparsity

    #     elif property_flag == VisualPropertyFlags.FIELD_AREA:
    #         return self.field_area

    #     elif property_flag == VisualPropertyFlags.COVERAGE:
    #         return self.coverage

    #     elif property_flag == VisualPropertyFlags.NUMEROSITY:
    #         return self.numerosity

    #     else:
    #         raise ValueError("f{property_flag} is a unknown visual feature")

    def to_dict(self) -> dict:
        """Dictionary with the visual properties"""
        rtn = [
            ("Hash", self._shapes.hash()),
            ("Numerosity", self.numerosity),
            (VisualPropertyFlags.AV_DOT_DIAMETER.label(), self.average_dot_diameter),
            (VisualPropertyFlags.AV_RECT_SIZE.label(),
             self.average_rectangle_size.tolist()),
            (VisualPropertyFlags.AV_PERIMETER.label(), self.average_perimeter),
            (VisualPropertyFlags.AV_SURFACE_AREA.label(), self.average_surface_area),
            (VisualPropertyFlags.TOTAL_PERIMETER.label(), self.total_perimeter),
            (VisualPropertyFlags.TOTAL_SURFACE_AREA.label(), self.total_surface_area),
            (VisualPropertyFlags.FIELD_AREA.label(), self.field_area),
            (VisualPropertyFlags.SPARSITY.label(), self.sparsity),
            (VisualPropertyFlags.COVERAGE.label(), self.coverage),
            (VisualPropertyFlags.LOG_SIZE.label(), self.log_size),
            (VisualPropertyFlags.LOG_SPACING.label(), self.log_spacing),
        ]
        return OrderedDict(rtn)
