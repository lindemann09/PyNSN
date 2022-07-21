__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod

import numpy as _np
from . import _colour
from .. import _lib
from .. import _arrays
from .._shapes import Dot, PictureFile

# helper for type checking and error raising error
def _check_object_array(obj):
    if not isinstance(obj, (_arrays.DotArray, _arrays.RectangleArray)):
        raise TypeError("DotArray or RectangleArray expected, but not {}".format(
            type(obj).__name__))


class ArrayDraw(metaclass=ABCMeta):
    """Generic array draw with abstract static methods

    To develop a plotter for other graphic system, inherit the abstract class
    and define you own drawing class (MyDraw())
        'get_squared_image', 'scale_image', 'draw_shape', 'draw_convex_hull'

    Image can be then generated via
    >>> MyDraw()().create_image(object_array=object_array, colours=colours)
    """

    @staticmethod
    @abstractmethod
    def get_squared_image(image_width, background_colour, **kwargs):
        """
        -------
        rtn : should return image
        """
        return

    @staticmethod
    @abstractmethod
    def scale_image(image, scaling_factor):
        pass

    @staticmethod
    @abstractmethod
    def draw_shape(image, shape, opacity, scaling_factor):
        """functions to draw object in the specific framework

        Returns
        -------
        image :  handler of plotter in the respective framework
                    (e.g. pillow image, axes (matplotlib) or svrdraw object)
        """
        pass

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
        pass


    def create_image(self, object_array, colours, antialiasing=None, **kwargs):
        """create image

        Parameters
        ----------
        object_array : the array
        colours : ImageColours
        antialiasing :   bool or number (scaling factor)
            Only useful for pixel graphics. If turn on, picture will be
            generated on a large pixel (cf. scaling factor) array and scaled
            down after generation


        Returns
        -------
        rtn : image
        """

        _check_object_array(object_array)
        if colours is None:
            colours = _colour.ImageColours()
        if not isinstance(colours, _colour.ImageColours):
            raise TypeError("Colours must be of type image.ImageColours")

        if isinstance(antialiasing, bool):
            if antialiasing:  # (not if 1)
                aaf = 2  # AA default
            else:
                aaf = 1
        else:
            try:
                aaf = int(antialiasing)
            except (ValueError, TypeError):
                aaf = 1

        # prepare the pil image, make target area if required
        image_width = int(_np.ceil(object_array.target_area_radius) * 2) * aaf
        img = self.get_squared_image(image_width=image_width,
                                     background_colour=colours.background.colour,
                                     **kwargs)

        if colours.target_area.colour is not None:
            obj = Dot(xy=(0, 0), diameter=image_width,
                           attribute=colours.target_area.colour)
            self.draw_shape(img, obj, opacity=1,
                            scaling_factor=1)

        if object_array.properties.numerosity > 0:

            # draw objects
            for obj in object_array.iter_objects():
                att = obj.get_attribute_object()
                if isinstance(att, PictureFile):
                    pass
                else:  # dot or rect: force colour, set default colour if no colour
                    obj.attribute = _colour.Colour(att, colours.default_object_colour)
                self.draw_shape(img, obj, opacity=colours.opacity_object,
                                scaling_factor=aaf)

            # draw convex hulls
            if colours.field_area_positions.colour is not None and \
                    object_array.properties.field_area_positions > 0:
                self.draw_convex_hull(img,
                                      points=object_array.properties.convex_hull_positions.xy,
                                      convex_hull_colour=colours.field_area_positions,
                                      opacity=colours.opacity_guides,
                                      scaling_factor=aaf)
            if colours.field_area.colour is not None and \
                    object_array.properties.field_area > 0:
                self.draw_convex_hull(img,
                                      points=object_array.properties.convex_hull.xy,
                                      convex_hull_colour=colours.field_area,
                                      opacity=colours.opacity_guides,
                                      scaling_factor=aaf)
            #  and center of mass
            if colours.center_of_field_area.colour is not None:
                obj = _lib.Dot(xy=object_array.get_center_of_field_area(),
                               diameter=10,
                               attribute=colours.center_of_field_area.colour)
                self.draw_shape(img, obj, opacity=colours.opacity_guides,
                                scaling_factor=aaf)
            if colours.center_of_mass.colour is not None:
                obj = _lib.Dot(xy=object_array.get_center_of_mass(),
                               diameter=10,
                               attribute=colours.center_of_mass.colour)
                self.draw_shape(img, obj, opacity=colours.opacity_guides,
                                scaling_factor=aaf)

        # rescale for antialiasing
        if aaf != 1:
            img = self.scale_image(img, scaling_factor=aaf)

        return img
