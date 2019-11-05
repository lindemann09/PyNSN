__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from PIL import Image, ImageDraw
import numpy as np
from .colour import ImageColours

def create(dot_array, colours, antialiasing=True,
           gabor_filter=None):
    # ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    antialiasing: Ture or integer

    gabor_filter: from PIL.ImageFilter
    default_dot_colour: if colour is undefined in dot_array
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

    if not isinstance(colours, ImageColours):
        raise ValueError("Colours must be a ImageColours instance")

    image_size = int(round(dot_array.target_array_radius * 2)) * aa
    img = Image.new("RGBA", (image_size, image_size),
                    color=colours.background.colour)

    tmp_colour = colours.target_area.colour
    if tmp_colour is not None:
        _draw_dot(img, xy=_convert_pos(np.zeros(2), image_size),
                  diameter=image_size,
                  colour=tmp_colour)

    # draw dots
    default_dot_colour = colours.default_dot_colour
    for xy, d, c in zip(_convert_pos(dot_array.rounded_xy * aa, image_size),
                        dot_array.diameters * aa,
                        dot_array.attributes.colours):
        if c.colour is None:
            c = default_dot_colour
        _draw_dot(img, xy=xy, diameter=d, colour=c.colour)  # todo draw pictures

    tmp_colour = colours.field_area.colour
    if tmp_colour is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull=_convert_pos(dot_array.convex_hull_positions * aa, image_size),
                          convex_hull_colour=tmp_colour)

    tmp_colour = colours.field_area_outer.colour
    if tmp_colour is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull=_convert_pos(dot_array.convex_hull_positions_full * aa, image_size),
                          convex_hull_colour=tmp_colour)

    tmp_colour = colours.center_of_mass.colour
    if tmp_colour is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_mass * aa, image_size),
                  diameter=10 * aa, colour=tmp_colour)

    tmp_colour = colours.center_of_outer_positions.colour
    if tmp_colour is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_outer_positions * aa, image_size),
                  diameter=10 * aa, colour=tmp_colour)

    if aa != 1:
        image_size = int(image_size / aa)
        img = img.resize((image_size, image_size), Image.LANCZOS)

    if gabor_filter is not None:
        try:
            img = img.filter(gabor_filter)
        except:
            raise RuntimeError("Can't apply gabor_filter {}".format(gabor_filter))

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
