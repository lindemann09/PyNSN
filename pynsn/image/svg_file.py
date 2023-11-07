__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import copy
from typing import Optional

import numpy as _np
import svgwrite as _svg

from .. import _stimulus, _shapes
from . import _base
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


class _SVGDraw(_base.AbstractArrayDraw):
    # scaling not used, because vector format is scale independent.

    @staticmethod
    def get_image(image_size, background_colour: str, **kwargs) -> _svg.Drawing:
        """"""
        size = (f"{image_size[0]}px", f"{image_size[1]}px")
        image = _svg.Drawing(size=size, filename=kwargs['filename'])
        if background_colour is not None:
            bkg_rect = _shapes.Rectangle(xy=(0, 0), size=image_size,
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
        if isinstance(shape, _shapes.Picture):
            raise RuntimeError("Pictures are not supported for SVG file.")

        shape = copy(shape)
        shape.xy = _base.cartesian2image_coordinates(shape.xy,
                                                     _np.array(svg_image_size(image))).tolist()
        col = shape.get_colour()

        if isinstance(shape, _shapes.Dot):
            image.add(image.circle(center=shape.xy,
                                   r=shape.diameter / 2,
                                   # stroke_width="0", stroke="black",
                                   fill=col.colour,
                                   opacity=opacity))
        elif isinstance(shape, _shapes.Rectangle):
            image.add(image.rect(insert=shape.left_bottom,
                                 size=shape.size,
                                 fill=col.colour,
                                 opacity=opacity))
        else:
            raise NotImplementedError(
                f"Shape {type(shape)} NOT YET IMPLEMENTED")

    @staticmethod
    def draw_convex_hull(image, points, convex_hull_colour, opacity,
                         scaling_factor):
        """"""

        points = _base.cartesian2image_coordinates(
            _np.asarray(points), _np.array(svg_image_size(image)))

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
