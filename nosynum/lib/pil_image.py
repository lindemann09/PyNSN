from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from PIL import Image, ImageDraw
import numpy as np


def create(dot_array,
           colour_area=None,
           colour_convex_hull_positions=None,
           colour_convex_hull_dots=None,
           colour_center_of_mass = None,
           colour_center_of_outer_positions=None,
           antialiasing=True,
           colour_background=(0, 0, 0)):
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image"""

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
    img = Image.new("RGBA", (image_size, image_size), color=colour_background)
    if colour_area is not None:
        _draw_dot(img, xy=_convert_pos(np.zeros(2), image_size),
                  diameter=image_size, colour=colour_area)

    for xy, d, c in zip(_convert_pos(dot_array.rounded_xy*aa, image_size),
                        dot_array.rounded_diameters*aa,
                        dot_array.colours):
        _draw_dot(img, xy=xy, diameter=d, colour=c)

    if colour_convex_hull_positions is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull = _convert_pos(dot_array.convex_hull_positions*aa, image_size),
                          convex_hull_colour=colour_convex_hull_positions)

    if colour_convex_hull_dots is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull = _convert_pos(dot_array.convex_hull_dots*aa, image_size),
                          convex_hull_colour=colour_convex_hull_dots)

    if colour_center_of_mass is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_mass*aa, image_size),
                  diameter=10*aa, colour=colour_center_of_mass)
    if colour_center_of_outer_positions is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_outer_positions*aa, image_size),
                  diameter=10*aa, colour=colour_center_of_outer_positions)

    if aa!=1:
        image_size = int(image_size / aa)
        img = img.resize((image_size, image_size), Image.LANCZOS)
    return img


def _draw_dot(img, xy, diameter, colour=None, picture=None):
    # draw a dot on an image

    if colour is None:
        colour = (255, 255, 255)
    else:
        colour = colour

    r = diameter // 2
    if picture is not None:
        pict = Image.open(picture, "r")
        img.paste(pict, (xy[0]-r, xy[1]-r))
    else:
        ImageDraw.Draw(img).ellipse((xy[0]-r, xy[1]-r, xy[0]+r, xy[1]+r), fill=tuple(colour))

def _convert_pos(xy, image_size):
    """convert dot pos to pil image coordinates"""
    return (xy * [1,-1]) + image_size//2

def _draw_convex_hull(img, convex_hull, convex_hull_colour):
    # plot convey hull
    hull = np.append(convex_hull, [convex_hull[0]], axis=0)
    last = None
    draw = ImageDraw.Draw(img)
    for p in hull:
        if last is not None:
            draw.line(np.append(last, p).tolist(),
                      width=2, fill=convex_hull_colour)
        last = p

