"""

"""

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import enum
from collections import OrderedDict
from typing import Any, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .._misc import key_value_format
from .._shapes import Dot, Ellipse, Picture, PolygonShape, Rectangle
from .._shapes import ellipse_geometry as ellipse_geo
from .shape_array import ShapeArray


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
        return (
            VisProp.LOG_SIZE,
            VisProp.TOTAL_SURFACE_AREA,
            VisProp.AV_SURFACE_AREA,
            VisProp.AV_PERIMETER,
            VisProp.TOTAL_PERIMETER,
        )

    @staticmethod
    def space_properties() -> tuple:
        return (VisProp.LOG_SPACING, VisProp.SPARSITY, VisProp.FIELD_AREA)

    def is_dependent_from(self, other: Any) -> bool:
        """returns true if both properties are not independent"""
        is_size_prop = self in VisProp.size_properties()
        is_space_prop = self in VisProp.space_properties()
        other_size_prop = other in VisProp.size_properties()
        other_space_prop = other in VisProp.space_properties()
        return (is_size_prop and other_size_prop) or (
            is_space_prop and other_space_prop
        )

    def __str__(self) -> str:
        if self == VisProp.NUMEROSITY:
            return "Numerosity"
        elif self == VisProp.LOG_SIZE:
            return "Log size"
        elif self == VisProp.TOTAL_SURFACE_AREA:
            return "Total surface area"
        elif self == VisProp.AV_SURFACE_AREA:
            return "Av. surface area"
        elif self == VisProp.AV_PERIMETER:
            return "Av. perimeter"
        elif self == VisProp.TOTAL_PERIMETER:
            return "Total perimeter"
        elif self == VisProp.LOG_SPACING:
            return "Log spacing"
        elif self == VisProp.SPARSITY:
            return "Sparsity"
        elif self == VisProp.FIELD_AREA:
            return "Field area"
        elif self == VisProp.COVERAGE:
            return "Coverage"
        else:
            return "???"

    def short_name(self) -> str:
        """Short names
        N=Numerosity
        TSA = Total surface area
        ASA = Av. surface area (also item surface area)
        AP = Av. perimeter
        TP=Total perimeter
        SP=Sparsity  (1/density)
        FA=Field area
        CO=Coverage
        """

        if self == VisProp.NUMEROSITY:
            return "N"
        elif self == VisProp.LOG_SIZE:
            return "logSize"
        elif self == VisProp.TOTAL_SURFACE_AREA:
            return "TSA"
        elif self == VisProp.AV_SURFACE_AREA:
            return "ASA"
        elif self == VisProp.AV_PERIMETER:
            return "AP"
        elif self == VisProp.TOTAL_PERIMETER:
            return "TP"
        elif self == VisProp.LOG_SPACING:
            return "logSpace"
        elif self == VisProp.SPARSITY:
            return "SP"
        elif self == VisProp.FIELD_AREA:
            return "FA"
        elif self == VisProp.COVERAGE:
            return "CO"
        else:
            return "???"


class ArrayProperties(object):
    """Non-Symbolic Number Stimulus"""

    def __init__(self, shape_array: ShapeArray) -> None:
        self._shape_arr = shape_array
        self._ch = None

    def to_text(self, short_format: bool = False) -> str:
        rtn = ""
        if not short_format:
            first = True
            for k, v in self.to_dict().items():
                if first and len(rtn) == 0:
                    rtn = "- "
                    first = False
                else:
                    rtn += " "
                rtn += key_value_format(k, v) + "\n "
        else:
            for k, v in self.to_dict(short_format=True).items():
                rtn += f"{k}: {v:.2f}, "
            rtn = rtn[:-2]
        return rtn.rstrip()

    def __repr__(self) -> str:
        return self.to_text(short_format=True)

    def __str__(self) -> str:
        return self.to_text()

    @property
    def areas(self) -> NDArray[np.float64]:
        """area of each object"""

        rtn = np.full(self._shape_arr.n_shapes, np.nan)
        # rects and polygons
        idx = np.append(
            self._shape_arr.get_ids(Rectangle.shape_type()),
            self._shape_arr.get_ids(Picture.shape_type())
        )
        if len(idx) > 0:
            rtn[idx] = self._shape_arr.sizes[idx, 0] * \
                self._shape_arr.sizes[idx, 1]

        # circular shapes area, Area = pi * r_x * r_y
        idx = np.append(
            self._shape_arr.get_ids(Dot.shape_type()),
            self._shape_arr.get_ids(Ellipse.shape_type()))
        if len(idx) > 0:
            r = self._shape_arr.sizes[idx, :] / 2
            rtn[idx] = np.pi * r[:, 0] * r[:, 1]

        # polygons area
        idx = self._shape_arr.get_ids(PolygonShape.shape_type())
        if len(idx) > 0:
            rtn[idx] = shapely.area(self._shape_arr.polygons[idx])
        return rtn

    @property
    def perimeter(self) -> NDArray[np.float64]:
        """Perimeter for each object"""

        rtn = np.full(self._shape_arr.n_shapes, np.nan)

        idx = np.concatenate((
            self._shape_arr.get_ids(Rectangle.shape_type()),
            self._shape_arr.get_ids(Picture.shape_type()),
            self._shape_arr.get_ids(PolygonShape.shape_type())
        ))
        if len(idx) > 0:
            rtn[idx] = shapely.length(self._shape_arr.polygons[idx])
        # dots perimeter
        idx = self._shape_arr.get_ids(Dot.shape_type())
        if len(idx) > 0:
            rtn[idx] = np.pi * self._shape_arr.sizes[idx, 0]
        # ellipse perimeter
        idx = self._shape_arr.get_ids(Ellipse.shape_type())
        if len(idx) > 0:
            rtn[idx] = ellipse_geo.perimeter(self._shape_arr.sizes[idx, :])

        return rtn

    @property
    def center_of_mass(self) -> NDArray:
        """center of mass of all shapes"""
        areas = self.areas
        weighted_sum = np.sum(self._shape_arr.xy *
                              np.atleast_2d(areas).T, axis=0)
        return weighted_sum / np.sum(areas)

    @property
    def numerosity(self) -> int:
        """number of shapes"""
        return self._shape_arr.n_shapes

    @property
    def total_surface_area(self) -> np.float64:
        return np.nansum(self.areas)

    @property
    def average_surface_area(self) -> np.float64:
        if self._shape_arr.n_shapes == 0:
            return np.float64(np.nan)
        return np.nanmean(self.areas)

    @property
    def total_perimeter(self) -> np.float64:
        return np.nansum(self.perimeter)

    @property
    def average_perimeter(self) -> np.float64:
        if self._shape_arr.n_shapes == 0:
            return np.float64(np.nan)
        return np.nanmean(self.perimeter)

    @property
    def coverage(self) -> np.float64:
        """percent coverage in the field area. It takes thus the object size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity
        """
        fa = self.field_area
        if fa == 0:
            return np.float64(np.nan)
        else:
            return self.total_surface_area / fa

    @property
    def log_size(self) -> np.float64:
        if self.numerosity == 0:
            return np.float64(np.nan)
        else:
            return np.log2(self.total_surface_area) + np.log2(self.average_surface_area)

    @property
    def log_spacing(self) -> np.float64:
        fa = self.field_area
        if fa == 0:
            return np.float64(np.nan)
        else:
            return np.log2(fa) + np.log2(self.sparsity)

    @property
    def sparsity(self) -> np.float64:
        if self.numerosity == 0:
            return np.float64(np.nan)
        else:
            return self.field_area / self.numerosity

    @property
    def field_area(self) -> np.float64:
        return np.float64(self._shape_arr.convex_hull.area)

    def get(self, prop: VisProp) -> Union[int, np.float64]:
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

    def to_dict(self, short_format: bool = False) -> dict:
        """Dictionary with the visual properties"""
        rtn = []
        if short_format:
            rtn.extend([(x.short_name(), self.get(x))
                        for x in list(VisProp)])  # type: ignore
        else:
            rtn.extend([(str(x), self.get(x))
                        for x in list(VisProp)])  # type: ignore
        return OrderedDict(rtn)
