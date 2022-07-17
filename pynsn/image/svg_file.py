__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
import svgwrite as _svg
from .._lib.geometry import cartesian2image_coordinates as _c2i_coord
from .. import _lib
from . import _array_draw


def create(object_array, colours, filename):
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


class _SVGDraw(_array_draw.ArrayDraw):
    # scaling not used, because vector format is scale independent.

    @staticmethod
    def get_squared_image(image_width, background_colour, **kwargs):
        """"""
        px = "{}px".format(image_width)
        return _svg.Drawing(size=(px, px), filename=kwargs['filename'])

    @staticmethod
    def scale_image(image, scaling_factor):
        """"""
        return image # not used

    @staticmethod
    def draw_shape(image, shape, opacity, scaling_factor):
        """"""
        assert isinstance(image, _svg.Drawing)
        width = int(image.attribs['width'][:-2]) # string "300px" --> 300
        shape.xy = _c2i_coord(_np.asarray(shape.xy), width).tolist()
        attr = shape.get_attribute_object()

        if isinstance(attr, _lib.PictureFile):
            raise RuntimeError("Pictures are not supported for SVG file.")

        if isinstance(shape, _lib.Dot):
            image.add(image.circle(center=shape.xy,
                                       r=shape.diameter / 2,
                                       # stroke_width="0", stroke="black",
                                       fill=attr.colour,
                                       opacity=opacity))
        elif isinstance(shape, _lib.Rectangle):
            image.add(image.rect(insert=(shape.left, shape.bottom),
                                     size=shape.size,
                                     fill=attr.colour,
                                     opacity=opacity))
        else:
            raise NotImplementedError("Shape {} NOT YET IMPLEMENTED".format(type(shape)))

    @staticmethod
    def draw_convex_hull(image, points, convex_hull_colour, opacity,
                         scaling_factor):
        """"""
        width = int(image.attribs['width'][:-2]) # string "300px" --> 300
        points = _c2i_coord(points, width)

        last = None
        for p in _np.append(points, [points[0]], axis=0):
            if last is not None:
                l = image.line(start=last, end=p).stroke(
                    width=1, color=convex_hull_colour.colour, opacity=opacity)
                image.add(l)
            last = p
