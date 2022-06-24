__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
import svgwrite as _svg
from . import _colour
from .._lib  import arrays as _arrays
from .._lib.geometry import cartesian2image_coordinates as _c2i_pos

def create(object_array, colours, filename="noname.svg"):
    assert isinstance(object_array, (_arrays.DotArray, _arrays.RectangleArray)) # FIXME implement for rect array
    if not isinstance(colours, _colour.ImageColours):
        raise TypeError("Colours must be of type pynsn.ImageColours")

    image_size = int(round(object_array.target_array_radius * 2))
    px = "{}px".format(image_size)
    svgdraw = _svg.Drawing(size = (px, px), filename=filename)

    if colours.target_area.colour is not None:
        svgdraw.add(svgdraw.circle(center=_c2i_pos(_np.zeros(2), image_size),
                                   r= image_size // 2,
                                   # stroke_width="0", stroke="black",
                                   fill=colours.target_area.colour))

    for xy, d, att in zip(_c2i_pos(object_array.xy, image_size),
                          object_array.diameters,
                          object_array.attributes):
        if att is None:
            c = colours.default_item_colour
        else:
            try:
                c = _colour.Colour(att)
            except TypeError:
                c = colours.default_item_colour

        svgdraw.add(svgdraw.circle(center=xy, r = d//2,
                                   #stroke_width="0", stroke="black",
                                    fill=c.colour))
    # FIXME TODO draw convex hulls
    return svgdraw
