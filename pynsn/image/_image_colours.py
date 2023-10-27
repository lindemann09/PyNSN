from collections import OrderedDict
from typing import Optional

from .. import _misc
from .._shapes import Colour
from ..types import ColourType


class ImageColours(object):
    COL_TARGET_AREA = "#FFF0D9"
    COL_DEFAULT_OBJECT = "darkgreen"
    OPACITY_OBJECT = 1
    OPACITY_GUIDES = 0.5

    def __init__(
        self,
        target_area: ColourType = None,
        field_area: ColourType = None,
        center_of_field_area: ColourType = None,
        center_of_mass: ColourType = None,
        background: ColourType = None,
        default_object_colour: ColourType = None,
        opacity_object: Optional[float] = None,
        opacity_guides: Optional[float] = None,
    ):
        if target_area is None:
            target_area = ImageColours.COL_TARGET_AREA
        if default_object_colour is None:
            default_object_colour = ImageColours.COL_DEFAULT_OBJECT
        if opacity_guides is None:
            opacity_guides = ImageColours.OPACITY_GUIDES
        if opacity_guides < 0 or opacity_guides > 1:
            raise ValueError(
                f"opacity_guides ({opacity_guides}) has to be between 0 and 1"
            )
        if opacity_object is None:
            opacity_object = ImageColours.OPACITY_OBJECT
        if opacity_object < 0 or opacity_object > 1:
            raise ValueError(
                f"opacity_object ({opacity_object}) has to be between 0 and 1"
            )

        self.target_area = Colour(target_area)
        self.field_area = Colour(field_area)
        self.center_of_field_area = Colour(center_of_field_area)
        self.center_of_mass = Colour(center_of_mass)
        self.background = Colour(background)
        self.default_object_colour = Colour(default_object_colour)
        self.opacity_object = opacity_object
        self.opacity_guides = opacity_guides

    def to_dict(self) -> dict:  # FIXME: is that used somewhere?
        return OrderedDict(
            {
                "total_area": self.target_area.value,
                "field_area": self.field_area.value,
                "center_of_field_area": self.center_of_field_area.value,
                "center_of_mass": self.center_of_mass.value,
                "background": self.background.value,
                "default_object": self.default_object_colour.value,
                "object_opacity": self.opacity_object,
                "info_opacity": self.opacity_guides,
            }
        )

    def __str__(self) -> str:
        return _misc.dict_to_text(self.to_dict())
