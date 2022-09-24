__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional, Union

import numpy as _np
from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw

from .._lib.geometry import cartesian2image_coordinates as _c2i_coord
from .. import _shapes
from . import _array_draw
from .. import _arrays
from ._image_colours import ImageColours


# TODO pillow supports no alpha/opacity


def create(object_array: _arrays.NSNStimulus,
           colours: Optional[ImageColours] = None,
           antialiasing: Union[bool, int] = True) -> _Image.Image:
    # ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    antialiasing: True or integer

    default_dot_colour: if colour is undefined in _lib
    """

    return _PILDraw().create_image(object_array=object_array,
                                   colours=colours,
                                   antialiasing=antialiasing)


class _PILDraw(_array_draw.ABCArrayDraw):

    @staticmethod
    def get_image(image_size, background_colour: str, **kwargs) -> _Image.Image:
        # filename not used for pil images
        return _Image.new("RGBA", image_size, color=background_colour)

    @staticmethod
    def scale_image(image, scaling_factor):
        try:
            resample = _Image.Resampling.LANCZOS  # type: ignore (old versions)
        except AttributeError:
            resample = _Image.LANCZOS

        im_size = (int(image.size[0] / scaling_factor),
                   int(image.size[1] / scaling_factor))
        return image.resize(im_size, resample=resample)

    @staticmethod
    def draw_shape(img, shape: _shapes.ShapeType,
                   opacity: float, scaling_factor: float):
        # FIXME opacity is ignored (not yet supported)
        # draw object
        shape.xy = _c2i_coord(shape.xy * scaling_factor, img.size).tolist()
        col = shape.get_colour()

        if isinstance(shape, _shapes.Dot):
            shape.diameter = shape.diameter * scaling_factor
            r = shape.diameter / 2
            _ImageDraw.Draw(img).ellipse((shape.x - r, shape.y - r,
                                          shape.x + r, shape.y + r),
                                         fill=col.colour)
        elif isinstance(shape, _shapes.Picture):
            tmp = _np.asarray(shape.size) * scaling_factor
            shape.size = tmp.tolist()
            # picture
            target_box = _np.round(shape.get_ltrb(), decimals=0)
            target_box[:, 1] = _np.flip(target_box[:, 1])  # reversed y axes
            pict = _Image.open(shape.filename, "r")
            if pict.size[0] != shape.size[0] or pict.size[1] != shape.size[1]:
                pict = pict.resize(shape.size.tolist(),
                                   resample=_Image.ANTIALIAS)

            tr_layer = _Image.new('RGBA', img.size, (0, 0, 0, 0))
            tr_layer.paste(pict, target_box)  # FIXME .flatten or .tolist()?
            res = _Image.alpha_composite(img, tr_layer)
            img.paste(res)

        elif isinstance(shape, _shapes.Rectangle):
            tmp = _np.asarray(shape.size) * scaling_factor
            shape.size = tmp.tolist()
            # rectangle shape
            _ImageDraw.Draw(img).rectangle((shape.left, shape.top,
                                            shape.right, shape.bottom),
                                           fill=col.colour)  # TODO decentral _shapes seems to be bit larger than with pyplot

        else:
            raise NotImplementedError(
                "Shape {} NOT YET IMPLEMENTED".format(type(shape)))

    @staticmethod
    def draw_convex_hull(img, points, convex_hull_colour,  opacity,
                         scaling_factor):
        # FIXME opacity is ignored (not yet supported)
        points = _c2i_coord(points * scaling_factor, img.size)
        last = None
        draw = _ImageDraw.Draw(img)
        for p in _np.append(points, [points[0]], axis=0):
            if last is not None:
                draw.line(_np.append(last, p).tolist(),
                          width=2,
                          fill=convex_hull_colour.colour)
            last = p
