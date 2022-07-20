__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw
import numpy as _np
from .._lib.geometry import cartesian2image_coordinates as _c2i_coord
from .. import _lib
from . import _array_draw


# TODO pillow supports no alpha/opacity


def create(object_array, colours=None, antialiasing=True):
    # ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    antialiasing: Ture or integer

    default_dot_colour: if colour is undefined in _lib
    """

    return _PILDraw().create_image(object_array=object_array,
                                   colours=colours,
                                   antialiasing=antialiasing)


class _PILDraw(_array_draw.ArrayDraw):

    @staticmethod
    def get_squared_image(image_width, background_colour, **kwargs):
        # filename not used for pil images
        return _Image.new("RGBA", (image_width, image_width), color=background_colour)

    @staticmethod
    def scale_image(image, scaling_factor):
        im_size = int(image.size[0] / scaling_factor)
        return image.resize((im_size, im_size), _Image.LANCZOS)

    @staticmethod
    def draw_shape(img, shape, opacity, scaling_factor):
        # FIXME opacity is ignored (not yet supported)
        # draw object
        shape.xy = _c2i_coord(_np.asarray(shape.xy) * scaling_factor, img.size[0]).tolist()
        attr = shape.get_attribute_object()

        if isinstance(shape, _lib.Dot):
            shape.diameter = shape.diameter * scaling_factor
            r = shape.diameter / 2
            _ImageDraw.Draw(img).ellipse((shape.x - r, shape.y - r,
                                          shape.x + r, shape.y + r),
                                         fill=attr.colour)
        elif isinstance(shape, _lib.Rectangle):
            tmp = _np.asarray(shape.size) * scaling_factor
            shape.size = tmp.tolist()
            if isinstance(attr, _lib.PictureFile):
                # picture
                shape_size = (round(shape.width), round(shape.height))
                target_box = (round(shape.left), round(shape.bottom),  # FIXME ceil?
                              round(shape.right), round(shape.top))  # reversed y axes
                pict = _Image.open(attr.filename, "r")
                if pict.size[0] != shape_size[0] or pict.size[1] != shape_size[1]:
                    pict = pict.resize(shape_size, resample=_Image.ANTIALIAS)

                tr_layer = _Image.new('RGBA', img.size, (0, 0, 0, 0))
                tr_layer.paste(pict, target_box)
                res = _Image.alpha_composite(img, tr_layer)
                img.paste(res)
            else:
                # rectangle shape
                _ImageDraw.Draw(img).rectangle((shape.left, shape.top,
                                                shape.right, shape.bottom),
                                               fill=attr.colour)  # FIXME decentral shapes seems to be bit larger than with pyplot

        else:
            raise NotImplementedError("Shape {} NOT YET IMPLEMENTED".format(type(shape)))

    @staticmethod
    def draw_convex_hull(img, points, convex_hull_colour,  opacity,
                         scaling_factor):
        # FIXME opacity is ignored (not yet supported)
        points = _c2i_coord(points * scaling_factor, img.size[0])
        last = None
        draw = _ImageDraw.Draw(img)
        for p in _np.append(points, [points[0]], axis=0):
            if last is not None:
                draw.line(_np.append(last, p).tolist(),
                          width=2,
                          fill=convex_hull_colour.colour)
            last = p
