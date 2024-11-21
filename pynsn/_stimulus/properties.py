"""

"""

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import enum
from collections import OrderedDict
from typing import Any, Dict, List, Tuple, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .._misc import key_value_format
from .._shapes import Dot, Ellipse, Picture, PolygonShape, Rectangle
from .._shapes import ellipse_geometry as ellipse_geo
from .shape_array import ShapeArray


class VisProp(enum.Enum):  # visual properties
    """Visual Properties

        `N` = Numerosity

        `TSA` = Total surface area

        `ASA` = Av. surface area or item surface area

        `AP` = Av. perimeter

        `TP` = Total perimeter

        `SP` = Sparsity  (=1/density)

        `FA` = Field area

        `CO` = Coverage

        `LOG_SIZE`

        `LOG_SPACING`
    """

    N = enum.auto()

    ASA = enum.auto()
    AP = enum.auto()

    TSA = enum.auto()
    TP = enum.auto()
    SP = enum.auto()
    FA = enum.auto()
    CO = enum.auto()

    LOG_SPACING = enum.auto()
    LOG_SIZE = enum.auto()

    @classmethod
    def space_properties(cls) -> Tuple[VisProp, VisProp, VisProp]:
        """tuple of all space properties"""
        return (cls.SP, cls.FA, cls.LOG_SPACING)

    @classmethod
    def size_properties(cls) -> Tuple[VisProp, VisProp, VisProp, VisProp, VisProp]:
        """tuple of all size properties"""
        return (cls.TSA, cls.ASA, cls.AP, cls.TP, cls.LOG_SIZE)

    @property
    def long_label(self) -> str:
        """long text representation of the visual property"""

        if self == VisProp.N:
            return "Numerosity"
        elif self == VisProp.LOG_SIZE:
            return "Log size"
        elif self == VisProp.TSA:
            return "Total surface area"
        elif self == VisProp.ASA:
            return "Av. item surface area"
        elif self == VisProp.AP:
            return "Av. item perimeter"
        elif self == VisProp.TP:
            return "Total perimeter"
        elif self == VisProp.LOG_SPACING:
            return "Log spacing"
        elif self == VisProp.SP:
            return "Sparsity"
        elif self == VisProp.FA:
            return "Field area"
        elif self == VisProp.CO:
            return "Coverage"
        else:
            return "???"

    def is_dependent_from(self, other: Any) -> bool:
        """returns true if both properties are not independent"""
        is_size_prop = self in VisProp.size_properties()
        is_space_prop = self in VisProp.space_properties()
        other_size_prop = other in VisProp.size_properties()
        other_space_prop = other in VisProp.space_properties()
        return (is_size_prop and other_size_prop) or (
            is_space_prop and other_space_prop)


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
            self._shape_arr.ids(Rectangle.shape_type()),
            self._shape_arr.ids(Picture.shape_type())
        )
        if len(idx) > 0:
            rtn[idx] = self._shape_arr.sizes[idx, 0] * \
                self._shape_arr.sizes[idx, 1]

        # circular shapes area, Area = pi * r_x * r_y
        idx = np.append(
            self._shape_arr.ids(Dot.shape_type()),
            self._shape_arr.ids(Ellipse.shape_type()))
        if len(idx) > 0:
            r = self._shape_arr.sizes[idx, :] / 2
            rtn[idx] = np.pi * r[:, 0] * r[:, 1]

        # polygons area
        idx = self._shape_arr.ids(PolygonShape.shape_type())
        if len(idx) > 0:
            rtn[idx] = shapely.area(self._shape_arr.polygons[idx])
        return rtn

    @property
    def perimeter(self) -> NDArray[np.float64]:
        """Perimeter for each object"""

        rtn = np.full(self._shape_arr.n_shapes, np.nan)

        idx = np.concatenate((
            self._shape_arr.ids(Rectangle.shape_type()),
            self._shape_arr.ids(Picture.shape_type()),
            self._shape_arr.ids(PolygonShape.shape_type())
        ))
        if len(idx) > 0:
            rtn[idx] = shapely.length(self._shape_arr.polygons[idx])
        # dots perimeter
        idx = self._shape_arr.ids(Dot.shape_type())
        if len(idx) > 0:
            rtn[idx] = np.pi * self._shape_arr.sizes[idx, 0]
        # ellipse perimeter
        idx = self._shape_arr.ids(Ellipse.shape_type())
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
        if prop == VisProp.AP:
            return self.average_perimeter

        elif prop == VisProp.TP:
            return self.total_perimeter

        elif prop == VisProp.ASA:
            return self.average_surface_area

        elif prop == VisProp.TSA:
            return self.total_surface_area

        elif prop == VisProp.LOG_SIZE:
            return self.log_size

        elif prop == VisProp.LOG_SPACING:
            return self.log_spacing

        elif prop == VisProp.SP:
            return self.sparsity

        elif prop == VisProp.FA:
            return self.field_area

        elif prop == VisProp.CO:
            return self.coverage

        elif prop == VisProp.N:
            return self.numerosity

        else:
            raise ValueError("f{property_flag} is a unknown visual feature")

    def to_dict(self, short_format: bool = False) -> dict:
        """Dictionary with the visual properties"""
        rtn = []
        if short_format:
            rtn.extend([(x.name, self.get(x))
                        for x in list(VisProp)])  # type: ignore
        else:
            rtn.extend([(x.long_label, self.get(x))
                        for x in list(VisProp)])  # type: ignore
        return OrderedDict(rtn)

# helper


VisPropList = List[str] | List[VisProp] | List[Union[str, VisProp]]


def ensure_vis_prop(prop: Union[str, VisProp]) -> VisProp:
    """helper fnc: creates a VisProp and raises if that is no possible"""

    if isinstance(prop, VisProp):
        return prop
    elif isinstance(prop, str):
        return VisProp[prop]
    else:
        raise ValueError(f"{prop} is not a visual feature.")


def ensure_vis_prop_list(prop: VisPropList) -> List[VisProp]:
    """helper fnc: creates a VisProp and raises if that is no possible"""

    rtn = []
    for p in prop:
        if isinstance(p, VisProp):
            rtn.append(p)
        elif isinstance(p, str):
            rtn.append(VisProp[p])
        else:
            raise ValueError(f"{p} is not a visual feature.")
    return rtn
