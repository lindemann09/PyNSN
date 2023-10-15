__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import copy
from typing import Optional

import numpy as _np
import svgwrite as _svg

from .. import _stimulus
from .._lib.geometry import cartesian2image_coordinates as _c2i_coord
from . import _array_draw
from ._image_colours import ImageColours


def create(filename: str,
           nsn_stimulus: _stimulus.NSNStimulus,
           colours: Optional[ImageColours] = None) -> _svg.Drawing:
    """SVG image/file, vector image format

    Parameters
    ----------
    nsn_stimulus
    colours
    filename

    Returns
    -------

    """
    return _SVGDraw().create_image(nsn_stimulus=nsn_stimulus, colours=colours,
                                   filename=filename)


class _SVGDraw(_array_draw.AbstractArrayDraw):
    # scaling not used, because vector format is scale independent.

    @staticmethod
    def get_image(image_size, background_colour: str, **kwargs) -> _svg.Drawing:
        """"""
        size = (f"{image_size[0]}px", f"{image_size[1]}px")
        image = _svg.Drawing(size=size, filename=kwargs['filename'])
        if background_colour is not None:
            bkg_rect = _stimulus.Rectangle(xy=(0, 0), size=image_size,
                                           attribute=background_colour)
            _SVGDraw.draw_shape(image=image, shape=bkg_rect, opacity=100,
                                scaling_factor=None)
        return image

    @staticmethod
    def scale_image(image, scaling_factor):
        """"""
        return image  # not used

    @staticmethod
    def draw_shape(image, shape, opacity, scaling_factor):
        """"""
        assert isinstance(image, _svg.Drawing)
        if isinstance(shape, _stimulus.Picture):
            raise RuntimeError("Pictures are not supported for SVG file.")

        shape = copy(shape)
        shape.xy = _c2i_coord(shape.xy,
                              _np.array(svg_image_size(image))).tolist()
        col = shape.get_colour()

        if isinstance(shape, _stimulus.Dot):
            image.add(image.circle(center=shape.xy,
                                   r=shape.diameter / 2,
                                   # stroke_width="0", stroke="black",
                                   fill=col.colour,
                                   opacity=opacity))
        elif isinstance(shape, _stimulus.Rectangle):
            image.add(image.rect(insert=(shape.left, shape.bottom),
                                 size=shape.size,
                                 fill=col.colour,
                                 opacity=opacity))
        else:
            raise NotImplementedError(
                "Shape {} NOT YET IMPLEMENTED".format(type(shape)))

    @staticmethod
    def draw_convex_hull(image, points, convex_hull_colour, opacity,
                         scaling_factor):
        """"""

        points = _c2i_coord(_np.asarray(points),
                            _np.array(svg_image_size(image)))

        last = None
        for p in _np.append(points, [points[0]], axis=0):
            if last is not None:
                l = image.line(start=last, end=p).stroke(
                    width=1, color=convex_hull_colour.colour, opacity=opacity)
                image.add(l)
            last = p


def svg_image_size(image):
    return (int(image.attribs['width'][:-2]),  # string "300px" --> 300
            int(image.attribs['height'][:-2]))
