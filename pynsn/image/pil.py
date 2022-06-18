__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw
import numpy as _np
from . import _colour
from .._lib.geometry import cartesian2image_coordinates as _c2i_coord


def create(dot_array, colours, antialiasing=True,
           gabor_filter=None):
    # ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    antialiasing: Ture or integer

    gabor_filter: from PIL.ImageFilter
    default_dot_colour: if colour is undefined in _nsn
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

    if not isinstance(colours, _colour.ImageColours):
        raise ValueError("Colours must be a ImageColours instance")

    image_size = int(round(dot_array.target_array_radius * 2)) * aa
    img = _Image.new("RGBA", (image_size, image_size),
                    color=colours.background.colour)

    if colours.target_area.colour is not None:
        _draw_dot(img, xy=_c2i_coord(_np.zeros(2), image_size),
                  diameter=image_size,
                  colour=colours.target_area.colour)

    dot_array = dot_array.copy()
    dot_array.round(decimals=0)

    # draw dots
    for xy, d, att in zip(_c2i_coord(dot_array.xy * aa, image_size),
                        dot_array.diameters * aa,
                          dot_array.attributes):
        if att is None:
            c = colours.default_dot_colour
        else:
            try:
                c = _colour.Colour(att)
            except:
                c = colours.default_dot_colour
        _draw_dot(img, xy=xy, diameter=d, colour=c.colour)

    if colours.field_area.colour is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull=_c2i_coord(
                              dot_array.features.convex_hull.xy * aa, image_size),
                          convex_hull_colour=colours.field_area.colour)

    if colours.field_area_outer.colour is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull=_c2i_coord(
                              dot_array.features.convex_hull.full_xy * aa,
                              image_size),
                          convex_hull_colour=colours.field_area_outer.colour)

    if colours.center_of_mass.colour is not None:
        _draw_dot(img, xy=_c2i_coord(dot_array.center_of_mass * aa, image_size),
                  diameter=10 * aa, colour=colours.center_of_mass.colour)

    if colours.center_of_outer_positions.colour is not None:
        _draw_dot(img, xy=_c2i_coord(dot_array.center_of_outer_positions * aa, image_size),
                  diameter=10 * aa, colour=colours.center_of_outer_positions.colour)

    if aa != 1:
        image_size = int(image_size / aa)
        img = img.resize((image_size, image_size), _Image.LANCZOS)

    if gabor_filter is not None:
        try:
            img = img.filter(gabor_filter)
        except:
            raise RuntimeError("Can't apply gabor_filter {}".format(gabor_filter))

    return img


def _draw_dot(img, xy, diameter, colour, picture=None):
    # draw a dot on an image

    r = diameter // 2
    if picture is not None:
        pict = _Image.open(picture, "r")
        img.paste(pict, (xy[0] - r, xy[1] - r))
    else:
        _ImageDraw.Draw(img).ellipse((xy[0] - r, xy[1] - r, xy[0] + r,
                                          xy[1] + r), fill=colour)


def _draw_convex_hull(img, convex_hull, convex_hull_colour):
    # plot convey hull

    hull = _np.append(convex_hull, [convex_hull[0]], axis=0)
    last = None
    draw = _ImageDraw.Draw(img)
    for p in hull:
        if last is not None:
            draw.line(_np.append(last, p).tolist(),
                      width=2,
                      fill=convex_hull_colour)
        last = p
