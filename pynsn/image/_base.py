__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from typing import Any, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._stimulus import NSNStimulus
from .._shapes import Dot, Picture, Rectangle, Ellipse

# helper for type checking and error raising error


def check_nsn_stimulus(obj):  # FIXME simpler (function not needed anymore)
    if not isinstance(obj, (NSNStimulus)):
        raise TypeError(
            "NSNStimulus expected, but not {}".format(type(obj).__name__))


def cartesian2image_coordinates(xy: ArrayLike,
                                image_size: ArrayLike) -> NDArray:
    """convert cartesian to image coordinates with (0,0) at top left and
    reversed y axis

    xy has to be a 2D array

    """
    return (np.asarray(xy) * [1, -1]) + np.asarray(image_size) / 2


class AbstractArrayDraw(metaclass=ABCMeta):
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
                    (e.g. pillow image, axes (matplotlib) or svgdraw object)
        """

    def create_image(
        self,
        nsn_stimulus: NSNStimulus,
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
        ta_shape = nsn_stimulus.target_area.shape
        colours = nsn_stimulus.colours
        if isinstance(ta_shape, Dot):
            target_area_shape = Dot(
                xy=(0, 0),
                diameter=np.ceil(ta_shape.diameter) * aaf,
                attribute=colours.target_area.value
            )
        elif isinstance(ta_shape, Ellipse):
            target_area_shape = Ellipse(
                xy=(0, 0),
                size=np.ceil(ta_shape.size) * aaf,
                attribute=colours.target_area.value
            )

        elif isinstance(ta_shape, Rectangle):
            target_area_shape = Rectangle(
                xy=(0, 0),
                size=np.ceil(ta_shape.size) * aaf,
                attribute=colours.target_area.value
            )
        else:
            raise NotImplementedError()  # should never happen

        image_size = (round(target_area_shape.width),
                      round(target_area_shape.height))
        img = self.get_image(
            image_size=image_size, background_colour=colours.background.value, **kwargs
        )

        if colours.target_area.value is not None:
            self.draw_shape(img, target_area_shape,
                            opacity=1, scaling_factor=1)

        if nsn_stimulus.properties.numerosity > 0:
            # draw objects
            for obj in nsn_stimulus.shapes:
                if obj.colour.value is None and not isinstance(obj, Picture):
                    # dot or rect: force colour, set default colour if no colour
                    obj.attribute = colours.object_default
                self.draw_shape(
                    img, obj, opacity=colours.opacity_object, scaling_factor=aaf)

            # draw convex hulls
            if colours.convex_hull.value is not None:
                coords = nsn_stimulus.convex_hull.coordinates
                if len(coords) > 1:
                    self.draw_convex_hull(img, points=coords,
                                          convex_hull_colour=colours.convex_hull,
                                          opacity=colours.opacity_guides,
                                          scaling_factor=aaf)
            #  and center of mass
            if colours.center_of_field_area.value is not None:
                obj = Dot(
                    xy=nsn_stimulus.convex_hull.centroid,
                    diameter=10,
                    attribute=colours.center_of_field_area,
                )
                self.draw_shape(
                    img, obj, opacity=colours.opacity_guides, scaling_factor=aaf
                )
            if colours.center_of_mass.value is not None:
                obj = Dot(
                    xy=nsn_stimulus.properties.center_of_mass,
                    diameter=10,
                    attribute=colours.center_of_mass,
                )
                self.draw_shape(
                    img, obj, opacity=colours.opacity_guides, scaling_factor=aaf
                )

        # rescale for antialiasing
        if aaf != 1:
            img = self.scale_image(img, scaling_factor=aaf)

        return img
