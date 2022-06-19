__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
import svgwrite as _svg
from . import _colour
from .._nsn.dot_array import DotArray as _DotArray
from .._lib.geometry import cartesian2image_coordinates as _c2i_pos

def create(dot_array, colours = _colour.ImageColours(), filename="noname.svg"):
    assert isinstance(dot_array, _DotArray)
    if not isinstance(colours, _colour.ImageColours):
        raise ValueError("Colours must be a pynsn.ImageColours instance")

    image_size = int(round(dot_array.target_array_radius * 2))
    px = "{}px".format(image_size)
    svgdraw = _svg.Drawing(size = (px, px), filename=filename)

    if colours.target_area.colour is not None:
        svgdraw.add(svgdraw.circle(center=_c2i_pos(_np.zeros(2), image_size),
                                   r= image_size // 2,
                                   # stroke_width="0", stroke="black",
                                   fill=colours.target_area.colour))

    dot_array = dot_array.copy()
    dot_array.round(decimals=1,int_type=float)
    for xy, d, att in zip(_c2i_pos(dot_array.xy, image_size),
                          dot_array.diameters,
                          dot_array.attributes):
        if att is None:
            c = colours.default_dot_colour
        else:
            try:
                c = _colour.Colour(att)
            except TypeError:
                c = colours.default_dot_colour

        svgdraw.add(svgdraw.circle(center=xy, r = d//2,
                                   #stroke_width="0", stroke="black",
                                    fill=c.colour))
    # FIXME TODO draw convex hulls
    return svgdraw
