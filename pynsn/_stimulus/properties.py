"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import enum
from collections import OrderedDict
from typing import Any, Union
import shapely
import numpy as np
from numpy.typing import NDArray

from .convex_hull import ConvexHull
from .shape_array import ShapeArray
from .shapes import Dot, Rectangle, Ellipse, Picture, PolygonShape
from .._lib.geometry import ellipse_perimeter


class VisProp(enum.Flag):  # visual properties
    NUMEROSITY = enum.auto()

    AV_SURFACE_AREA = enum.auto()
    AV_PERIMETER = enum.auto()

    TOTAL_SURFACE_AREA = enum.auto()
    TOTAL_PERIMETER = enum.auto()
    SPARSITY = enum.auto()
    FIELD_AREA = enum.auto()
    COVERAGE = enum.auto()

    LOG_SPACING = enum.auto()
    LOG_SIZE = enum.auto()

    @staticmethod
    def size_properties() -> tuple:
        return (VisProp.LOG_SIZE,
                VisProp.TOTAL_SURFACE_AREA,
                VisProp.AV_SURFACE_AREA,
                VisProp.AV_PERIMETER,
                VisProp.TOTAL_PERIMETER)

    @staticmethod
    def space_properties() -> tuple:
        return (VisProp.LOG_SPACING,
                VisProp.SPARSITY,
                VisProp.FIELD_AREA)

    def is_dependent_from(self, other: Any) -> bool:
        """returns true if both properties are not independent"""
        is_size_prop = self in VisProp.size_properties()
        is_space_prop = self in VisProp.space_properties()
        other_size_prop = other in VisProp.size_properties()
        other_space_prop = other in VisProp.space_properties()
        return (is_size_prop and other_size_prop) or \
            (is_space_prop and other_space_prop)

    def label(self) -> str:
        labels = {
            VisProp.NUMEROSITY: "Numerosity",
            VisProp.LOG_SIZE: "Log size",
            VisProp.TOTAL_SURFACE_AREA: "Total surface area",
            VisProp.AV_SURFACE_AREA: "Av. surface area",
            VisProp.AV_PERIMETER: "Av. perimeter",
            VisProp.TOTAL_PERIMETER: "Total perimeter",
            VisProp.LOG_SPACING: "Log spacing",
            VisProp.SPARSITY: "Sparsity",
            VisProp.FIELD_AREA: "Field area",
            VisProp.COVERAGE: "Coverage"}
        return labels[self]


class ArrayProperties(object):
    """Non-Symbolic Number Stimulus"""

    def __init__(self, shape_array: ShapeArray) -> None:
        self._shapes = shape_array
        self._ch = None

    @property
    def areas(self) -> NDArray:
        """area of each object"""

        rtn = np.full(self._shapes.n_objects, np.nan)
        # rects and polygons
        idx = np.append(self._shapes.ids[Rectangle.ID],
                        self._shapes.ids[Picture.ID])
        if len(idx) > 0:
            rtn[idx] = self._shapes.sizes[idx, 0] * self._shapes.sizes[idx, 1]
        # circular shapes area
        # Area = pi * r_x * r_y
        idx = np.append(self._shapes.ids[Dot.ID],
                        self._shapes.ids[Ellipse.ID])
        if len(idx) > 0:
            r = self._shapes.sizes[idx, :] / 2
            rtn[idx] = np.pi * r[:, 0] * r[:, 1]
        # polygons area
        idx = self._shapes.ids[PolygonShape.ID]
        if len(idx) > 0:
            rtn[idx] = shapely.area(self._shapes.polygons[idx])
        return rtn

    @property
    def perimeter(self) -> NDArray:
        """Perimeter for each dot"""

        rtn = np.full(self._shapes.n_objects, np.nan)
        # rects and polygons
        idx = np.append(self._shapes.ids[Rectangle.ID],
                        self._shapes.ids[Picture.ID],
                        self._shapes.ids[PolygonShape.ID])
        if len(idx) > 0:
            rtn[idx] = shapely.length(self._shapes.polygons[idx])
        # dots perimeter
        idx = self._shapes.ids[Dot.ID]
        if len(idx) > 0:
            rtn[idx] = np.pi * self._shapes.sizes[idx, 0]
        # ellipse perimeter
        idx = self._shapes.ids[Ellipse.ID]
        if len(idx) > 0:
            rtn[idx] = ellipse_perimeter(self._shapes.sizes[idx, :])

        return rtn

    @property
    def center_of_mass(self) -> NDArray:
        """center of mass of all objects"""
        areas = self.areas
        weighted_sum = np.sum(self._shapes.xy * np.atleast_2d(areas).T, axis=0)
        return weighted_sum / np.sum(areas)

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
    def total_surface_area(self) -> float:
        return np.nansum(self.areas)

    @property
    def average_surface_area(self) -> float:
        if self._shapes.n_objects == 0:
            return 0.0
        return np.nanmean(self.areas)

    @property
    def total_perimeter(self) -> float:
        return np.nansum(self.perimeter)

    @property
    def average_perimeter(self) -> float:
        if self._shapes.n_objects == 0:
            return 0.0
        return np.nanmean(self.perimeter)

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

    def get(self, prop: VisProp) -> Any:
        """returns a visual property"""
        if prop == VisProp.AV_PERIMETER:
            return self.average_perimeter

        elif prop == VisProp.TOTAL_PERIMETER:
            return self.total_perimeter

        elif prop == VisProp.AV_SURFACE_AREA:
            return self.average_surface_area

        elif prop == VisProp.TOTAL_SURFACE_AREA:
            return self.total_surface_area

        elif prop == VisProp.LOG_SIZE:
            return self.log_size

        elif prop == VisProp.LOG_SPACING:
            return self.log_spacing

        elif prop == VisProp.SPARSITY:
            return self.sparsity

        elif prop == VisProp.FIELD_AREA:
            return self.field_area

        elif prop == VisProp.COVERAGE:
            return self.coverage

        elif prop == VisProp.NUMEROSITY:
            return self.numerosity

        else:
            raise ValueError("f{property_flag} is a unknown visual feature")

    def to_dict(self) -> dict:
        """Dictionary with the visual properties"""
        rtn = []
        rtn.extend([(x.label(), self.get(x))
                   for x in list(VisProp)])  # type: ignore
        return OrderedDict(rtn)
