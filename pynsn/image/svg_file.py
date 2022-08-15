__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional

import numpy as _np
import svgwrite as _svg

from .. import _arrays, _shapes
from .._lib.np_coordinates import cartesian2image_coordinates as _c2i_coord
from . import _array_draw, _colour


def create(filename: str,
           object_array: _arrays.ObjectArrayType,
           colours: Optional[_colour.ImageColours] = None) -> _svg.Drawing:
    """SVG image/file, vector image format

    Parameters
    ----------
    object_array
    colours
    filename

    Returns
    -------

    """
    return _SVGDraw().create_image(object_array=object_array, colours=colours,
                                   filename=filename)


class _SVGDraw(_array_draw.ABCArrayDraw):
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
        shape.xy = _c2i_coord(_np.asarray(shape.xy),
                              _np.array(svg_image_size(image))).tolist()
        attr = shape.get_attribute_object()
        if isinstance(attr, _shapes.PictureFile):
            raise RuntimeError("Pictures are not supported for SVG file.")

        if isinstance(shape, _shapes.Dot):
            image.add(image.circle(center=shape.xy,
                                   r=shape.diameter / 2,
                                   # stroke_width="0", stroke="black",
                                   fill=attr.colour,
                                   opacity=opacity))
        elif isinstance(shape, _shapes.Rectangle):
            size = (float(shape.width), float(shape.height))
            image.add(image.rect(insert=(shape.left, shape.bottom),
                                 size=size,
                                 fill=attr.colour,
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
