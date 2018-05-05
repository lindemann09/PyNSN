from __future__ import print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from PIL import Image, ImageDraw
import numpy as np
from .._lib.colour import Colour
from .._lib.misc import PYTHON3


class PILImageGenerator(object):

    def __init__(self,
                 colour_target_area=None,
                 colour_field_area=None,
                 colour_field_area_outer=None,
                 colour_center_of_mass=None,
                 colour_center_of_outer_positions=None,
                 antialiasing=None,
                 colour_background=None,
                 default_dot_colour=Colour("lightgreen")  # used if no color specified in dot array
                 ):

        self.colour_target_area = Colour(colour_target_area)
        self.colour_field_area = Colour(colour_field_area)
        self.colour_field_area_outer = Colour(colour_field_area_outer)
        self.colour_center_of_mass = Colour(colour_center_of_mass)
        self.colour_center_of_outer_positions = Colour(colour_center_of_outer_positions)
        self.colour_background = Colour(colour_background)
        self.default_dot_colour = Colour(default_dot_colour)

        if isinstance(antialiasing, bool):
            if antialiasing:  # (not if 1)
                self.antialiasing = 2  # AA default
            else:
                self.antialiasing = 1
        else:
            try:
                self.antialiasing = int(antialiasing)
            except:
                self.antialiasing = 1

    def as_dict(self):
        return {"colour_total_area": self.colour_target_area.colour,
                "colour_field_area": self.colour_field_area.colour,
                "colour_field_area_outer": self.colour_field_area_outer.colour,
                "colour_center_of_mass": self.colour_center_of_mass.colour,
                "colour_center_of_outer_positions": self.colour_center_of_outer_positions.colour,
                "antialiasing": self.antialiasing,
                "colour_background": self.colour_background.colour,
                "default_dot_colour": self.default_dot_colour.colour}

    def make(self, dot_array):  # todo using *args and ImageParameter
        """use PIL colours (see PIL.ImageColor.colormap)

        returns pil image

        default_dot_colour: if coulor is undefined in dot_array
        """

        aa = self.antialiasing
        image_size = int(round(dot_array.target_array_radius * 2)) * aa
        img = Image.new("RGBA", (image_size, image_size),
                        color=self.colour_background.colour)

        tmp_colour = self.colour_target_area.colour
        if tmp_colour is not None:
            _draw_dot(img, xy=_convert_pos(np.zeros(2), image_size),
                      diameter=image_size,
                      colour=tmp_colour)

        # draw dots
        default_dot_colour = self.default_dot_colour
        for xy, d, c in zip(_convert_pos(dot_array.rounded_xy * aa, image_size),
                            dot_array.rounded_diameters * aa,
                            dot_array.attributes.colours):
            if c.colour is None:
                c = default_dot_colour
            _draw_dot(img, xy=xy, diameter=d, colour=c.colour)  # todo draw pictures

        tmp_colour = self.colour_field_area.colour
        if tmp_colour is not None:
            # plot convey hull
            _draw_convex_hull(img=img,
                              convex_hull=_convert_pos(dot_array.convex_hull_positions * aa, image_size),
                              convex_hull_colour=tmp_colour)

        tmp_colour = self.colour_field_area_outer.colour
        if tmp_colour is not None:
            # plot convey hull
            _draw_convex_hull(img=img,
                              convex_hull=_convert_pos(dot_array.full_convex_hull_positions * aa, image_size),
                              convex_hull_colour=tmp_colour)

        tmp_colour = self.colour_center_of_mass.colour
        if tmp_colour is not None:
            _draw_dot(img, xy=_convert_pos(dot_array.center_of_mass * aa, image_size),
                      diameter=10 * aa, colour=tmp_colour)

        tmp_colour = self.colour_center_of_outer_positions.colour
        if tmp_colour is not None:
            _draw_dot(img, xy=_convert_pos(dot_array.center_of_outer_positions * aa, image_size),
                      diameter=10 * aa, colour=tmp_colour)

        if aa != 1:
            image_size = int(image_size / aa)
            img = img.resize((image_size, image_size), Image.LANCZOS)

        self._image = img
        return img


def _convert_pos(xy, image_size):
    """convert dot pos to pil image coordinates"""
    return (xy * [1, -1]) + image_size // 2


def _draw_dot(img, xy, diameter, colour, picture=None):
    # draw a dot on an image

    r = diameter // 2
    if picture is not None:
        pict = Image.open(picture, "r")
        img.paste(pict, (xy[0] - r, xy[1] - r))
    else:
        ImageDraw.Draw(img).ellipse((xy[0] - r, xy[1] - r, xy[0] + r, xy[1] + r), fill=colour)


def _draw_convex_hull(img, convex_hull, convex_hull_colour):
    # plot convey hull

    hull = np.append(convex_hull, [convex_hull[0]], axis=0)
    last = None
    draw = ImageDraw.Draw(img)
    for p in hull:
        if last is not None:
            draw.line(np.append(last, p).tolist(),
                      width=2,
                      fill=convex_hull_colour)
        last = p
