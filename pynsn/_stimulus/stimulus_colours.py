from collections import OrderedDict

from .. import _misc, defaults
from .._shapes.colour import Colour, ColourLike


class StimulusColours(object):
    def __init__(
        self,
        target_area: ColourLike = defaults.COLOUR_TARGET_AREA,
        object_default: ColourLike = defaults.COLOUR_OBJECT,
        background: ColourLike = None,
        convex_hull: ColourLike = None,
        center_of_field_area: ColourLike = None,
        center_of_mass: ColourLike = None,
        opacity_object: float = defaults.OPACITY_OBJECT,
        opacity_guides: float = defaults.OPACITY_GUIDES
    ):
        self._target_area = Colour(target_area)
        self._convex_hull = Colour(convex_hull)
        self._object_default = Colour(object_default)
        self._center_of_field_area = Colour(center_of_field_area)
        self._center_of_mass = Colour(center_of_mass)
        self._background = Colour(background)
        self._opacity_object = 0.0
        self._opacity_guides = 0.0
        # call setter
        self.opacity_object = opacity_object
        self.opacity_guides = opacity_guides

    def todict(self) -> dict:
        return OrderedDict(
            {
                "total_area": self.target_area.value,
                "background": self.background.value,
                "default_object": self.object_default.value,
                "convex_hull": self.convex_hull.value,
                "center_of_field_area": self.center_of_field_area.value,
                "center_of_mass": self.center_of_mass.value,
                "object_opacity": self.opacity_object,
                "info_opacity": self.opacity_guides,
            }
        )

    def __str__(self) -> str:
        return _misc.dict_to_text(self.todict())

    @property
    def target_area(self) -> Colour:
        """ """
        return self._target_area

    @target_area.setter
    def target_area(self, val: ColourLike):
        self._target_area = Colour(val)

    @property
    def convex_hull(self) -> Colour:
        """ """
        return self._convex_hull

    @convex_hull.setter
    def convex_hull(self, val: ColourLike):
        self._convex_hull = Colour(val)

    @property
    def center_of_field_area(self) -> Colour:
        """ """
        return self._center_of_field_area

    @center_of_field_area.setter
    def center_of_field_area(self, val: ColourLike):
        self._center_of_field_area = Colour(val)

    @property
    def center_of_mass(self) -> Colour:
        """ """
        return self._center_of_mass

    @center_of_mass.setter
    def center_of_mass(self, val: ColourLike):
        self._center_of_mass = Colour(val)

    @property
    def background(self) -> Colour:
        """ """
        return self._background

    @background.setter
    def background(self, val: ColourLike):
        self._background = Colour(val)

    @property
    def object_default(self) -> Colour:
        """ """
        return self._object_default

    @object_default.setter
    def object_default(self, val: ColourLike):
        self._object_default = Colour(val)

    @property
    def opacity_object(self) -> float:
        """ """
        return self._opacity_object

    @opacity_object.setter
    def opacity_object(self, val: float):
        if val < 0 or val > 1:
            raise ValueError("opacity_object has to be between 0 and 1")
        self._opacity_object = val

    @property
    def opacity_guides(self) -> float:
        """ """
        return self._opacity_guides

    @opacity_guides.setter
    def opacity_guides(self, val: float):
        if val < 0 or val > 1:
            raise ValueError("opacity_guides has to be between 0 and 1")
        self._opacity_guides = val
