
import numpy as _np
from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw
from pynsn import Colour, _shapes
from pynsn._lib.np_coordinates import cartesian2image_coordinates as _c2i_coord


class Line(object):
    """to draw lines in the test pictures"""

    def __init__(self, xy_a, xy_b, colour):
        self.xy_a = xy_a
        self.xy_b = xy_b
        self.colour = Colour(colour).colour


def _draw_shape(img, shape, scaling_factor=1):

    if isinstance(shape, Line):
        xy_a = _c2i_coord(_np.asarray(shape.xy_a), img.size).tolist()
        xy_b = _c2i_coord(_np.asarray(shape.xy_b), img.size).tolist()
        xy_a.extend(xy_b)
        _ImageDraw.Draw(img).line(xy_a, fill=shape.colour, width=2)
        return

    shape.xy = _c2i_coord(_np.asarray(shape.xy) *
                          scaling_factor, img.size).tolist()
    attr = shape.get_attribute_object()
    if isinstance(shape, _shapes.Dot):
        shape.diameter = shape.diameter * scaling_factor
        r = shape.diameter / 2
        _ImageDraw.Draw(img).ellipse((shape.x - r, shape.y - r,
                                      shape.x + r, shape.y + r),
                                     fill=attr.colour)
    elif isinstance(shape, _shapes.Rectangle):
        tmp = _np.asarray(shape.size) * scaling_factor
        shape.size = tmp.tolist()
        if isinstance(attr, _shapes.PictureFile):
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
                                           fill=attr.colour)  # FIXME decentral _shapes seems to be bit larger than with pyplot

    else:
        raise NotImplementedError(
            "Shape {} NOT YET IMPLEMENTED".format(type(shape)))


def shapes_test_picture(shapes,  size=(500, 500),
                        filename="shapes_test.png",
                        background_colour="#888888"):
    """makes the shape test picture"""

    img = _Image.new("RGBA", size, color=background_colour)
    for s in shapes:
        _draw_shape(img, s)
    img.save(filename)
    print(f"saved {filename}")
