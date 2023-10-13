__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from typing import Any, Optional

import numpy as _np

from .._lib.colour import Colour
from .._shapes import Dot, Rectangle, Picture
from .._stimulus import NSNStimulus
from ._image_colours import ImageColours

# helper for type checking and error raising error


def check_nsn_stimulus(obj):  # FIXME simpler (function not needed anymore)
    if not isinstance(obj, (NSNStimulus)):
        raise TypeError("NSNStimulus expected, but not {}".format(type(obj).__name__))


class ABCArrayDraw(metaclass=ABCMeta):
    """Generic array draw with abstract static methods

    To develop a plotter for other graphic system, inherit the abstract class
    and define you own drawing class (MyDraw())
        'get_image', 'scale_image', 'draw_shape', 'draw_convex_hull'

    Image can be then generated via
    >>> MyDraw()().create_image(nsn_stimulus=nsn_stimulus, colours=colours)
    """

    @staticmethod
    @abstractmethod
    def get_image(image_size, background_colour, **kwargs) -> Any:
        """
        -------
        rtn : should return image
        """
        return

    @staticmethod
    @abstractmethod
    def scale_image(image, scaling_factor):
        """ """

    @staticmethod
    @abstractmethod
    def draw_shape(image, shape, opacity, scaling_factor):
        """functions to draw object in the specific framework

        Returns
        -------
        image :  handler of plotter in the respective framework
                    (e.g. pillow image, axes (matplotlib) or svrdraw object)
        """

    @staticmethod
    @abstractmethod
    def draw_convex_hull(image, points, convex_hull_colour, opacity, scaling_factor):
        """functions to draw object in the specific framework

        Parameters
        ----------
        opacity
        scaling_factor
        convex_hull_colour
        points
        image :  handler of plotter in the respective framework
                    (e.g. pillow image, axes (matplotlib) or svrdraw object)
        """

    def create_image(
        self,
        nsn_stimulus: NSNStimulus,
        colours: Optional[ImageColours],
        antialiasing: Optional[float] = None,
        **kwargs
    ) -> Any:
        """create image

        Parameters
        ----------
        nsn_stimulus : the array
        colours : ImageColours
        antialiasing :   bool or number (scaling factor)
            Only useful for pixel graphics. If turn on, picture will be
            generated on a large pixel (cf. scaling factor) array and scaled
            down after generation


        Returns
        -------
        rtn : image
        """

        check_nsn_stimulus(nsn_stimulus)
        if colours is None:
            colours = ImageColours()
        if not isinstance(colours, ImageColours):
            raise TypeError("Colours must be of type image.ImageColours")

        if isinstance(antialiasing, bool):
            if antialiasing:  # (not if 1)
                aaf = 2  # AA default
            else:
                aaf = 1
        else:
            try:
                aaf = int(antialiasing)  # type: ignore
            except (ValueError, TypeError):
                aaf = 1

        # prepare the image, make target area if required
        if isinstance(nsn_stimulus.target_area_shape, Dot):
            tmp = int(_np.ceil(nsn_stimulus.target_area_shape.diameter) * aaf)
            target_area_shape = Dot(
                xy=(0, 0), diameter=tmp, attribute=colours.target_area.value
            )
            image_size = _np.ones(2) * tmp

        elif isinstance(nsn_stimulus.target_area_shape, Rectangle):
            tmp = _np.int16(_np.ceil(nsn_stimulus.target_area_shape.size) * aaf)
            target_area_shape = Rectangle(
                xy=(0, 0), size=tmp, attribute=colours.target_area.value
            )
            image_size = target_area_shape.size
        else:
            raise NotImplementedError()  # should never happen

        image_size = (round(image_size[0]), round(image_size[1]))
        img = self.get_image(
            image_size=image_size, background_colour=colours.background.value, **kwargs
        )

        if colours.target_area.value is not None:
            self.draw_shape(img, target_area_shape, opacity=1, scaling_factor=1)

        if nsn_stimulus.properties.numerosity > 0:
            # draw objects
            for obj in nsn_stimulus.objects.iter():
                att = obj.get_colour()
                if att.colour is None and not isinstance(obj, Picture):
                    # dot or rect: force colour, set default colour if no colour
                    obj.attribute = Colour(None, colours.default_object_colour.value)
                self.draw_shape(
                    img, obj, opacity=colours.opacity_object, scaling_factor=aaf
                )

            # draw convex hulls
            if (
                colours.field_area_positions.value is not None
                and nsn_stimulus.properties.field_area_positions > 0
            ):
                self.draw_convex_hull(
                    img,
                    points=nsn_stimulus.properties.convex_hull_positions.xy,
                    convex_hull_colour=colours.field_area_positions,
                    opacity=colours.opacity_guides,
                    scaling_factor=aaf,
                )
            if (
                colours.field_area.value is not None
                and nsn_stimulus.properties.field_area > 0
            ):
                self.draw_convex_hull(
                    img,
                    points=nsn_stimulus.properties.convex_hull.xy,
                    convex_hull_colour=colours.field_area,
                    opacity=colours.opacity_guides,
                    scaling_factor=aaf,
                )
            #  and center of mass
            if colours.center_of_field_area.value is not None:
                obj = Dot(
                    xy=nsn_stimulus.properties.convex_hull.center,
                    diameter=10,
                    attribute=colours.center_of_field_area.value,
                )
                self.draw_shape(
                    img, obj, opacity=colours.opacity_guides, scaling_factor=aaf
                )
            if colours.center_of_mass.value is not None:
                obj = Dot(
                    xy=nsn_stimulus.objects.center_of_mass,
                    diameter=10,
                    attribute=colours.center_of_mass.value,
                )
                self.draw_shape(
                    img, obj, opacity=colours.opacity_guides, scaling_factor=aaf
                )

        # rescale for antialiasing
        if aaf != 1:
            img = self.scale_image(img, scaling_factor=aaf)

        return img
