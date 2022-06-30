__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw
import numpy as _np
from . import _colour
from .._lib.geometry import cartesian2image_coordinates as _c2i_coord
from .._lib import shapes as _shape
from .._lib  import arrays as _arrays

from ._colour import ImageColours # make available

#TODO pillow supports no alpha

def create(object_array, colours, antialiasing=True, _gabor_filter=None):
    # ImageParameter
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image

    antialiasing: Ture or integer

    gabor_filter: from PIL.ImageFilter
    default_dot_colour: if colour is undefined in _lib
    """

    _arrays._check_object_array(object_array)
    assert isinstance(colours, _colour.ImageColours)

    if isinstance(antialiasing, bool):
        if antialiasing:  # (not if 1)
            aaf = 2  # AA default
        else:
            aaf = 1
    else:
        try:
            aaf = int(antialiasing)
        except:
            aaf = 1

    # prepare the pil image, make target area if required
    image_size = int(_np.ceil(object_array.target_area_radius) * 2) * aaf
    img = _Image.new("RGBA", (image_size, image_size),
                     color=colours.background.colour)

    if colours.target_area.colour is not None:
        obj = _shape.Dot(xy=_c2i_coord(_np.zeros(2), image_size),
                         diameter=image_size,
                         attribute=colours.target_area.colour)
        _draw_shape(img, obj)

    if object_array.features.numerosity > 0:
        image_coord = _c2i_coord(object_array.xy * aaf, image_size)
        if isinstance(object_array, _arrays.DotArray):
            # draw dots
            for xy, d, att in zip(image_coord, object_array.diameters * aaf,
                                  object_array.attributes):
                obj = _shape.Dot(xy=xy, diameter=d)
                obj.attribute = _colour.make_colour(att,
                                                    colours.default_object_colour)
                _draw_shape(img, obj)

        elif isinstance(object_array, _arrays.RectangleArray):
            # draw rectangle
            for xy, size, att in zip(image_coord,
                                     object_array.sizes * aaf,
                                     object_array.attributes):
                obj = _shape.Rectangle(xy=xy, size=size)
                obj.attribute = _colour.make_colour(att,
                                                    colours.default_object_colour)
                _draw_shape(img, obj)

        # draw convex hulls
        if colours.field_area_positions.colour is not None:
            _draw_convex_hull(img=img,
                              points=_c2i_coord(
                                  object_array.features.convex_hull_positions.xy * aaf, image_size),
                              convex_hull_colour=colours.field_area_positions.colour)
        if colours.field_area.colour is not None:
            _draw_convex_hull(img=img,
                              points=_c2i_coord(
                                  object_array.features.convex_hull.xy * aaf,
                                  image_size),
                              convex_hull_colour=colours.field_area.colour)
        #  and center of mass
        if colours.center_of_positions.colour is not None:
            obj = _shape.Dot(xy=_c2i_coord(object_array.center_of_mass() * aaf, image_size),
                             diameter=10 * aaf,
                             attribute=colours.center_of_positions.colour)
            _draw_shape(img, obj)
        if colours.center_of_mass.colour is not None:
            obj = _shape.Dot(xy=_c2i_coord(object_array.center_of_positions() * aaf, image_size),
                             diameter=10 * aaf,
                             attribute=colours.center_of_mass.colour)
            _draw_shape(img, obj)

    # rescale for antialising
    if aaf != 1:
        image_size = int(image_size / aaf)
        img = img.resize((image_size, image_size), _Image.LANCZOS)

    # TODO gabor needed?
    if _gabor_filter is not None:
        try:
            img = img.filter(_gabor_filter)
        except:
            raise RuntimeError("Can't apply gabor_filter {}".format(_gabor_filter))

    return img


def _draw_shape(img, shape):
    # draw object
    assert isinstance(shape, (_shape.Dot, _shape.Rectangle))

    colour = _colour.Colour(shape.attribute)
    if isinstance(shape, _shape.Dot):
        r = shape.diameter / 2
        _ImageDraw.Draw(img).ellipse((shape.x - r, shape.y - r,
                                      shape.x + r, shape.y + r),
                                     fill=colour.colour)

    elif isinstance(shape, _shape.Rectangle):
        _ImageDraw.Draw(img).rectangle((shape.left, shape.top,
                                        shape.right, shape.bottom),
                                       fill=colour.colour) # FIXME decentral shapes a bit large than with pyplot

    else:
        raise NotImplementedError("Shape {} NOT YET IMPLEMENTED".format(type(shape)))

    # TODO pictures in attributes
    #if picture is not None:
    #    pict = _Image.open(picture, "r")
    #    img.paste(pict, (xy[0] - r, xy[1] - r))

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

