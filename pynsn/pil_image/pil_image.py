from __future__ import print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from PIL import Image, ImageDraw
import numpy as np
from .._lib.colour import Colour


def create(dot_array,
           colour_area=None,
           colour_convex_hull_positions=None,
           colour_convex_hull_dots=None,
           colour_center_of_mass = None,
           colour_center_of_outer_positions=None,
           antialiasing=True,
           colour_background=(0, 0, 0),
           default_dot_colour="lightgreen"):
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    default_dot_colour: if coulor is undefined in dot_array
    """

    if isinstance(antialiasing, bool):
        if antialiasing: # (not if 1)
            aa = 2 # AA default
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

    #draw dots
    default_dot_colour = Colour(default_dot_colour)
    for xy, d, c in zip(_convert_pos(dot_array.rounded_xy*aa, image_size),
                        dot_array.rounded_diameters*aa,
                        dot_array.features.colours):
        if c.colour is None:
            c = default_dot_colour
        _draw_dot(img, xy=xy, diameter=d, colour=c.colour) # todo draw pictures

    tmp_colour = Colour(colour_convex_hull_positions).colour
    if tmp_colour is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull = _convert_pos(dot_array.convex_hull_positions*aa, image_size),
                          convex_hull_colour = tmp_colour)

    tmp_colour = Colour(colour_convex_hull_dots).colour
    if tmp_colour is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull = _convert_pos(dot_array.convex_hull * aa, image_size),
                          convex_hull_colour = tmp_colour)

    tmp_colour = Colour(colour_center_of_mass).colour
    if tmp_colour is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_mass*aa, image_size),
                  diameter=10*aa, colour = tmp_colour)

    tmp_colour = Colour(colour_center_of_outer_positions).colour
    if tmp_colour is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_outer_positions*aa, image_size),
                  diameter=10*aa, colour = tmp_colour)

    if aa!=1:
        image_size = int(image_size / aa)
        img = img.resize((image_size, image_size), Image.LANCZOS)

    return img

def _convert_pos(xy, image_size):
    """convert dot pos to pil image coordinates"""
    return (xy * [1,-1]) + image_size//2


def _draw_dot(img, xy, diameter, colour, picture=None):
    # draw a dot on an image

    r = diameter // 2
    if picture is not None:
        pict = Image.open(picture, "r")
        img.paste(pict, (xy[0]-r, xy[1]-r))
    else:
        ImageDraw.Draw(img).ellipse((xy[0]-r, xy[1]-r, xy[0]+r, xy[1]+r), fill=colour)

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

