from __future__ import print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from collections import namedtuple
from PIL import Image, ImageDraw
import numpy as np
from .._lib.colour import Colour
from .._lib.generator import DotArrayGenerator


class ImageParameter(object):

    def __init__(self, colour_area=None,
                 colour_convex_hull_positions=None,
                 colour_convex_hull_dots=None,
                 colour_center_of_mass=None,
                 colour_center_of_outer_positions=None,
                 antialiasing=None,
                 colour_background=None,
                 default_dot_colour="lightgreen"):
        self.colour_area = colour_area
        self.colour_convex_hull_positions = colour_convex_hull_positions
        self.colour_convex_hull_dots = colour_convex_hull_dots
        self.colour_center_of_mass = colour_center_of_mass
        self.colour_center_of_outer_positions = colour_center_of_outer_positions
        self.antialiasing = antialiasing
        self.colour_background = colour_background
        self.default_dot_colour = default_dot_colour # used if no color specified in dot array

#number
# da_generator: max_array_radius dot_colour dot_diameter_mean dot_diameter_range dot_diameter_std minimum_gap

RandomDotImageParameter = namedtuple("DotPixmapParameter",  # all paramter of da generator and pil image
                                     "number max_array_radius dot_colour dot_diameter_mean " +
                                     "dot_diameter_range dot_diameter_std minimum_gap colour_area " +
                                     "colour_convex_hull_positions colour_convex_hull_dots colour_center_of_mass " +
                                     "colour_center_of_outer_positions antialiasing colour_background")
RandomDotImageParameter.__new__.__defaults__ = (None,) * len(RandomDotImageParameter._fields) # TODO use *args method

def create(dot_array,
           colour_area=None,
           colour_convex_hull_positions=None,
           colour_convex_hull_dots=None,
           colour_center_of_mass=None,
           colour_center_of_outer_positions=None,
           antialiasing=True,
           colour_background=(0, 0, 0),
           default_dot_colour="lightgreen"): #todo using *args and ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    default_dot_colour: if coulor is undefined in dot_array
    """

    if isinstance(antialiasing, bool):
        if antialiasing:  # (not if 1)
            aa = 2  # AA default
        else:
            aa = 1
    else:
        try:
            aa = int(antialiasing)
        except:
            aa = 1

    image_size = int(round(dot_array.max_array_radius * 2)) * aa
    img = Image.new("RGBA", (image_size, image_size),
                    color=Colour(colour_background).colour)

    tmp_colour = Colour(colour_area).colour
    if tmp_colour is not None:
        _draw_dot(img, xy=_convert_pos(np.zeros(2), image_size),
                  diameter=image_size,
                  colour=tmp_colour)

    # draw dots
    default_dot_colour = Colour(default_dot_colour)
    for xy, d, c in zip(_convert_pos(dot_array.rounded_xy * aa, image_size),
                        dot_array.rounded_diameters * aa,
                        dot_array.features.colours):
        if c.colour is None:
            c = default_dot_colour
        _draw_dot(img, xy=xy, diameter=d, colour=c.colour)  # todo draw pictures

    tmp_colour = Colour(colour_convex_hull_positions).colour
    if tmp_colour is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull=_convert_pos(dot_array.convex_hull_positions * aa, image_size),
                          convex_hull_colour=tmp_colour)

    tmp_colour = Colour(colour_convex_hull_dots).colour
    if tmp_colour is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull=_convert_pos(dot_array.convex_hull * aa, image_size),
                          convex_hull_colour=tmp_colour)

    tmp_colour = Colour(colour_center_of_mass).colour
    if tmp_colour is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_mass * aa, image_size),
                  diameter=10 * aa, colour=tmp_colour)

    tmp_colour = Colour(colour_center_of_outer_positions).colour
    if tmp_colour is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_outer_positions * aa, image_size),
                  diameter=10 * aa, colour=tmp_colour)

    if aa != 1:
        image_size = int(image_size / aa)
        img = img.resize((image_size, image_size), Image.LANCZOS)

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


# class PILImageGenerator(DotArrayGenerator, ImageParameter): TODO


def generate_random_dot_array_image(para, logger=None):
    """
    Generate randam Dor Array from RandomDotImageParameter
    para: RandomDotImageParameter
    logger has to be a GeneratorLogger

    returns image and dot_array
    """

    if not isinstance(para, RandomDotImageParameter):
        raise TypeError("para has to be RandomDotImageParameter, but not {}".format(
                                type(para).__name__))

    generator = DotArrayGenerator(
        max_array_radius=para.max_array_radius,
        dot_diameter_mean=para.dot_diameter_mean,
        dot_diameter_range=para.dot_diameter_range,
        dot_diameter_std=para.dot_diameter_std,
        dot_colour=para.dot_colour,
        minimum_gap=para.minimum_gap        )

    dot_array = generator.make(n_dots=para.number, logger=logger)
    image = create(dot_array,
                   colour_area=para.colour_area,
                   colour_convex_hull_positions=para.colour_convex_hull_positions,
                   colour_convex_hull_dots=para.colour_convex_hull_dots,
                   colour_center_of_mass=para.colour_center_of_mass,
                   colour_center_of_outer_positions=para.colour_center_of_outer_positions,
                   antialiasing=para.antialiasing,
                   colour_background=para.colour_background)

    return image, dot_array
