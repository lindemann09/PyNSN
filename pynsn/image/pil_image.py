"""
"""
__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import typing as _tp

import numpy as _np
from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw

from .._lib .geometry import cartesian2image_coordinates as _c2i_coord
from . import _array_draw
from .. import _stimulus
from ._image_colours import ImageColours

# TODO pillow supports no alpha/opacity

RESAMPLING = _Image.Resampling.LANCZOS


def create(
    nsn_stimulus: _stimulus.NSNStimulus,
    colours: _tp.Optional[ImageColours] = None,
    antialiasing: _tp.Union[bool, int] = True,
) -> _Image.Image:
    # ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    antialiasing: True or integer

    default_dot_colour: if colour is undefined in _lib
    """

    return _PILDraw().create_image(
        nsn_stimulus=nsn_stimulus, colours=colours, antialiasing=antialiasing
    )


class _PILDraw(_array_draw.AbstractArrayDraw):
    @staticmethod
    def get_image(image_size, background_colour: str, **kwargs) -> _Image.Image:
        # filename not used for pil images
        return _Image.new("RGBA", image_size, color=background_colour)

    @staticmethod
    def scale_image(image, scaling_factor):
        im_size = (
            int(image.size[0] / scaling_factor),
            int(image.size[1] / scaling_factor),
        )
        return image.resize(im_size, resample=RESAMPLING)

    @staticmethod
    def draw_shape(
        image, shape: _stimulus.ShapeType, opacity: float, scaling_factor: float
    ):
        # FIXME opacity is ignored (not yet supported)
        # draw object
        shape.copy()
        shape.xy = _c2i_coord(_np.asarray(shape.xy) *
                              scaling_factor, image.size)

        if isinstance(shape, _stimulus.Dot):
            r = (shape.diameter * scaling_factor) / 2
            x, y = shape.xy
            _ImageDraw.Draw(image).ellipse(
                (x - r, y - r, x + r, y + r), fill=shape.colour.value
            )

        elif isinstance(shape, _stimulus.Picture):
            shape.size = _np.asarray(shape.size) * scaling_factor
            # picture
            target_box = _np.array([shape.left, shape.top,
                                   shape.right, shape.bottom])
            target_box[:, 1] = _np.flip(target_box[:, 1])  # reversed y axes
            pict = _Image.open(shape.filename, "r")
            if pict.size[0] != shape.size[0] or pict.size[1] != shape.size[1]:
                pict = pict.resize((int(shape.size[0]), int(shape.size[1])),
                                   resample=RESAMPLING)

            tr_layer = _Image.new("RGBA", image.size, (0, 0, 0, 0))
            tr_layer.paste(pict, target_box.flatten().tolist())
            res = _Image.alpha_composite(image, tr_layer)
            image.paste(res)

        elif isinstance(shape, _stimulus.Rectangle):
            shape.size = _np.asarray(shape.size) * scaling_factor
            # rectangle shape
            _ImageDraw.Draw(image).rectangle((shape.left, shape.bottom,
                                              shape.right, shape.top),
                                             fill=shape.colour.value)  # TODO decentral _shapes seems to be bit larger than with pyplot

        else:
            raise NotImplementedError(
                f"Shape {type(shape)} NOT YET IMPLEMENTED")

    @staticmethod
    def draw_convex_hull(image, points, convex_hull_colour, opacity, scaling_factor):
        # FIXME opacity is ignored (not yet supported)
        points = _c2i_coord(points * scaling_factor, image.size)
        last = None
        draw = _ImageDraw.Draw(image)
        for p in _np.append(points, [points[0]], axis=0):
            if last is not None:
                draw.line(_np.append(last, p).tolist(),
                          width=2, fill=convex_hull_colour.value)
            last = p
