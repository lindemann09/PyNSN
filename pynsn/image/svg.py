__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
import svgwrite as _svg
from ..lib.colour import ImageColours as _ImageColours


def create(dot_array, colours = _ImageColours(), filename="noname.svg"):
    if not isinstance(colours, _ImageColours):
        raise ValueError("Colours must be a ImageColours instance")

    image_size = int(round(dot_array.target_array_radius * 2))
    px = "{}px".format(image_size)
    svgdraw = _svg.Drawing(size = (px, px), filename=filename)

    if colours.target_area.colour is not None:
        svgdraw.add(svgdraw.circle(center=_convert_pos(_np.zeros(2), image_size),
                                   r= image_size // 2,
                                   # stroke_width="0", stroke="black",
                                   fill=colours.target_area.colour))

    dot_array = dot_array.copy()
    dot_array.round(decimals=1,int_type=float)
    for xy, d, c in zip(_convert_pos(dot_array.xy, image_size),
                        dot_array.diameters,
                        dot_array.get_colours()):
        if c.colour is None:
            c = colours.default_dot_colour
        svgdraw.add(svgdraw.circle(center=xy, r = d//2,
                                   #stroke_width="0", stroke="black",
                                    fill=c.colour))

    # FIXME TODO draw convex hulls


    return svgdraw


def _convert_pos(xy, image_size): # TODO coordinate system
    """convert dot pos to svg coordinates"""
    return (xy * [1, -1]) + image_size // 2