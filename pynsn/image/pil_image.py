__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw
import numpy as _np
from . import _colour
from .._lib.geometry import cartesian2image_coordinates as _c2i_coord
from .. import shapes as _shapes
from .. import arrays as _arrays
from .._lib import misc as _misc


from ._colour import ImageColours # make available

#FIXME pillow supports no alpha

def create(object_array, colours=None, antialiasing=True, _gabor_filter=None):
    # ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    antialiasing: Ture or integer

    gabor_filter: from PIL.ImageFilter
    default_dot_colour: if colour is undefined in _lib
    """

    _arrays._check_object_array(object_array)
    if colours is None:
        colours = _colour.ImageColours()
    if not isinstance(colours, _colour.ImageColours):
        raise TypeError("Colours must be of type image.ImageColours")

    if isinstance(antialiasing, bool):
        if antialiasing:  # (not if 1)
            aaf = 2  # AA default
        else:
            aaf = 1
    else:
        try:
            aaf = int(antialiasing)
        except ValueError:
            aaf = 1

    # prepare the pil image, make target area if required
    image_size = int(_np.ceil(object_array.target_area_radius) * 2) * aaf
    img = _Image.new("RGBA", (image_size, image_size),
                     color=colours.background.colour)

    if colours.target_area.colour is not None:
        obj = _shapes.Dot(xy=_c2i_coord(_np.zeros(2), image_size),
                          diameter=image_size,
                          attribute=colours.target_area.colour)
        _draw_shape(img, obj)

    if object_array.properties.numerosity > 0:
        image_coord = _c2i_coord(object_array.xy * aaf, image_size)
        if isinstance(object_array, _arrays.DotArray):
            # draw dots
            for xy, d, att in zip(image_coord, object_array.diameters * aaf,
                                  object_array.attributes):
                obj = _shapes.Dot(xy=xy, diameter=d)
                obj.attribute = _colour.Colour(att,
                                               colours.default_object_colour)
                _draw_shape(img, obj)

        elif isinstance(object_array, _arrays.RectangleArray):
            # draw rectangle
            for xy, size, att in zip(image_coord,
                                     object_array.sizes * aaf,
                                     object_array.attributes):
                obj = _shapes.Rectangle(xy=xy, size=size,
                                        attribute=att)
                if obj.get_picture() is None:
                    # rect
                    obj.attribute = _colour.Colour(att,
                                colours.default_object_colour)
                _draw_shape(img, obj)

        # draw convex hulls
        if colours.field_area_positions.colour is not None and \
                object_array.properties.field_area_positions > 0:
            _draw_convex_hull(img=img,
                              points=_c2i_coord(
                                  object_array.properties.convex_hull_positions.xy * aaf, image_size),
                              convex_hull_colour=colours.field_area_positions.colour)
        if colours.field_area.colour is not None and \
                object_array.properties.field_area > 0:
            _draw_convex_hull(img=img,
                              points=_c2i_coord(
                                  object_array.properties.convex_hull.xy * aaf,
                                  image_size),
                              convex_hull_colour=colours.field_area.colour)
        #  and center of mass
        if colours.center_of_field_area.colour is not None:
            obj = _shapes.Dot(xy=_c2i_coord(object_array.center_of_field_area() * aaf, image_size),
                              diameter=10 * aaf,
                              attribute=colours.center_of_field_area.colour)
            _draw_shape(img, obj)
        if colours.center_of_mass.colour is not None:
            obj = _shapes.Dot(xy=_c2i_coord(object_array.center_of_mass() * aaf, image_size),
                              diameter=10 * aaf,
                              attribute=colours.center_of_mass.colour)
            _draw_shape(img, obj)

    # rescale for antialiasing
    if aaf != 1:
        image_size = int(image_size / aaf)
        img = img.resize((image_size, image_size), _Image.LANCZOS)

    # FIXME gabor needed?
    if _gabor_filter is not None:
        try:
            img = img.filter(_gabor_filter)
        except:
            raise RuntimeError("Can't apply gabor_filter {}".format(_gabor_filter))

    return img


def _draw_shape(img, shape):
    # draw object
    colour = _colour.Colour(shape.attribute)
    if isinstance(shape, _shapes.Dot):
        r = shape.diameter / 2
        _ImageDraw.Draw(img).ellipse((shape.x - r, shape.y - r,
                                      shape.x + r, shape.y + r),
                                     fill=colour.colour)
    elif isinstance(shape, _shapes.Rectangle):
        if shape.get_picture() is not None:
            # picture
            shape_size = (round(shape.width), round(shape.height))
            target_box = (round(shape.left), round(shape.bottom),
                          round(shape.right), round(shape.top)) # reversed y axes
            pict = _Image.open(shape.get_picture(), "r")
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
                                       fill=colour.colour) # FIXME decentral shapes seems to be bit larger than with pyplot

    else:
        raise NotImplementedError("Shape {} NOT YET IMPLEMENTED".format(type(shape)))


def _draw_convex_hull(img, points, convex_hull_colour):
    # plot convey hull

    last = None
    draw = _ImageDraw.Draw(img)
    for p in _np.append(points, [points[0]], axis=0):
        if last is not None:
            draw.line(_np.append(last, p).tolist(),
                      width=2,
                      fill=convex_hull_colour)
        last = p

