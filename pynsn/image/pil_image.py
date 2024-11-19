"""
"""
__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import math as _math
import typing as _tp
from copy import deepcopy

import numpy as _np
from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw

from .. import NSNStimulus as _NSNStimulus
from .. import NSNStimulusPair as _NSNStimulusPair
from .. import _shapes
from ..defaults import VERT_ECCENTRICITY_STIMPAIR
from . import _base

# TODO pillow supports no alpha/opacity

RESAMPLING = _Image.Resampling.LANCZOS


def create(
    nsn_stimulus: _NSNStimulus,
    antialiasing: _tp.Union[bool, int] = True,
) -> _Image.Image:
    # ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    antialiasing: True or integer

    """

    return _PILDraw().create_image(
        nsn_stimulus=nsn_stimulus, antialiasing=antialiasing
    )


def create_stim_pair(stim_pair: _NSNStimulusPair,
                     postion_a: _tp.Optional[_tp.Tuple[int, int]] = None,
                     postion_b: _tp.Optional[_tp.Tuple[int, int]] = None,
                     swap_positions: bool = False,
                     padding: int = 10,
                     antialiasing: _tp.Union[bool, int] = True,
                     background_image: _tp.Optional[_Image.Image] = None) -> _Image.Image:
    """returns a pil image with two NSNStimuli, one left and one right

    Note
    ----
    see create
    """
    if postion_a is None:
        postion_a = (-VERT_ECCENTRICITY_STIMPAIR, 0)
    if postion_b is None:
        postion_b = (+VERT_ECCENTRICITY_STIMPAIR, 0)

    if swap_positions:
        tmp = postion_a
        postion_a = postion_b
        postion_b = tmp

    im = (create(stim_pair.stim_a, antialiasing),
          create(stim_pair.stim_b, antialiasing))

    max_width = max(im[0].size[0], im[1].size[0])
    max_height = max(im[0].size[1], im[1].size[1])
    max_abs_x = max(abs(postion_a[0]), abs(postion_b[0]))
    max_abs_y = max(abs(postion_a[1]), abs(postion_b[1]))
    # height, width, center
    bkg_size2 = (max_abs_x + _math.ceil(max_width/2) + padding,
                 max_abs_y + _math.ceil(max_height/2) + padding)

    if isinstance(background_image, _Image.Image):
        bkg = background_image
    else:
        # (0, 0, 0, 0) is fully transparent
        bkg = _Image.new("RGBA",
                         (bkg_size2[0]*2, bkg_size2[1]*2), (0, 0, 0, 0))

    bkg.paste(im[0], (postion_a[0] + bkg_size2[0] - im[0].size[0]//2,
                      postion_a[1] + bkg_size2[1] - im[0].size[1]//2))
    bkg.paste(im[1], (postion_b[0] + bkg_size2[0] - im[1].size[0]//2,
                      postion_b[1] + bkg_size2[1] - im[1].size[1]//2))
    return bkg


class _PILDraw(_base.AbstractArrayDraw):
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
        image, shape: _shapes.AbstractShape, opacity: float, scaling_factor: float
    ):
        # FIXME opacity is ignored (not yet supported)
        # draw object
        shape = deepcopy(shape)
        shape.xy = _base.cartesian2image_coordinates(
            _np.asarray(shape.xy) * scaling_factor, image.size)
        shape.scale(scaling_factor)

        if isinstance(shape, (_shapes.Ellipse, _shapes.Dot)):
            rx, ry = shape.size / 2
            x, y = shape.xy
            _ImageDraw.Draw(image).ellipse(
                (x - rx, y - ry, x + rx, y + ry), fill=shape.colour.value
            )

        elif isinstance(shape, _shapes.Picture):
            upper_left = _np.flip(shape.left_top).tolist()
            pict = _Image.open(shape.path, "r")
            if pict.size[0] != shape.size[0] or pict.size[1] != shape.size[1]:
                pict = pict.resize((int(shape.size[0]), int(shape.size[1])),
                                   resample=RESAMPLING)

            tr_layer = _Image.new("RGBA", image.size, (0, 0, 0, 0))
            tr_layer.paste(pict, upper_left)
            res = _Image.alpha_composite(image, tr_layer)
            image.paste(res)

        elif isinstance(shape, _shapes.Rectangle):
            # rectangle shape TODO decentral _shapes seems to be bit larger than with pyplot
            _ImageDraw.Draw(image).rectangle(tuple(shape.box),  # type: ignore
                                             fill=shape.colour.value)
        else:
            raise NotImplementedError(
                f"Shape {type(shape)} NOT YET IMPLEMENTED")

    @staticmethod
    def draw_convex_hull(image, points, convex_hull_colour, opacity, scaling_factor):
        # FIXME opacity is ignored (not yet supported)
        points = _base.cartesian2image_coordinates(
            points * scaling_factor, image.size)
        last = None
        draw = _ImageDraw.Draw(image)
        for p in _np.append(points, [points[0]], axis=0):
            if last is not None:
                draw.line(_np.append(last, p).tolist(),
                          width=2, fill=convex_hull_colour.value)
            last = p
