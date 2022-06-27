__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
import svgwrite as _svg
from . import _colour
from .._lib  import arrays as _arrays
from .._lib.geometry import cartesian2image_coordinates as _c2i_coord
from .._lib import shape as _shape

def create(object_array, colours, filename="noname.svg"):
    assert isinstance(object_array, (_arrays.DotArray, _arrays.RectangleArray)) # FIXME implement for rect array
    if not isinstance(colours, _colour.ImageColours):
        raise TypeError("Colours must be of type pynsn.ImageColours")

    image_size = _np.ceil(object_array.target_array_radius) * 2
    px = "{}px".format(image_size)
    svgdraw = _svg.Drawing(size = (px, px), filename=filename)

    if colours.target_area.colour is not None:
        svgdraw.add(svgdraw.circle(center=_c2i_coord(_np.zeros(2), image_size),
                                   r= object_array.target_array_radius,
                                   # stroke_width="0", stroke="black",
                                   fill=colours.target_area.colour))

    if object_array.features.numerosity > 0:
        image_coord = _c2i_coord(object_array.xy, image_size)
        if isinstance(object_array, _arrays.DotArray):
            # draw dots
            for xy, d, att in zip(image_coord, object_array.diameters,
                                  object_array.attributes):
                obj = _shape.Dot(xy=xy, diameter=d)
                obj.attribute = _colour.make_colour(att,
                                                    colours.default_object_colour)
                _draw_shape(svgdraw, obj, opacity=colours.object_opacity)

        elif isinstance(object_array, _arrays.RectangleArray):
            # draw rectangle
            for xy, size, att in zip(image_coord,
                                     object_array.sizes,
                                     object_array.attributes):
                obj = _shape.Rectangle(xy=xy, size=size)
                obj.attribute = _colour.make_colour(att,
                                            colours.default_object_colour)
                _draw_shape(svgdraw, obj, opacity=colours.object_opacity)

        # draw convex hulls
        if colours.field_area_position.colour is not None:
            _draw_convex_hull(svgdraw=svgdraw,
                      points=_c2i_coord(
                          object_array.features.convex_hull.position_xy, image_size),
                      convex_hull_colour=colours.field_area_position.colour,
                      opacity=colours.info_shapes_opacity)
        if colours.field_area_outer.colour is not None:
            _draw_convex_hull(svgdraw=svgdraw,
                      points=_c2i_coord(
                          object_array.features.convex_hull.outer_xy,
                          image_size),
                      convex_hull_colour=colours.field_area_outer.colour,
                      opacity=colours.info_shapes_opacity)
        #  and center of mass
        if colours.center_of_mass.colour is not None:
            obj = _shape.Dot(xy=_c2i_coord(object_array.center_of_mass(), image_size),
                             diameter=10,
                             attribute=colours.center_of_mass.colour)
            _draw_shape(svgdraw, obj, opacity=colours.info_shapes_opacity)
        if colours.center_of_outer_positions.colour is not None:
            obj = _shape.Dot(xy=_c2i_coord(object_array.center_of_outer_positions, image_size),
                             diameter=10,
                             attribute=colours.center_of_outer_positions.colour)
            _draw_shape(svgdraw, obj, opacity=colours.info_shapes_opacity)

    return svgdraw


def _draw_shape(svgdraw, shape, opacity=1):
    # draw object
    assert isinstance(shape, (_shape.Dot, _shape.Rectangle))

    colour = _colour.Colour(shape.attribute)

    if isinstance(shape, _shape.Dot):
        r = shape.diameter / 2
        svgdraw.add(svgdraw.circle(center=shape.xy,
                                   r=shape.diameter/2,
                                   # stroke_width="0", stroke="black",
                                   fill=colour.colour,
                                   opacity=opacity))

    elif isinstance(shape, _shape.Rectangle):
        svgdraw.add(svgdraw.rect(insert=(shape.left, shape.bottom),
                                 size=shape.size,
                                 fill=colour.colour,
                                 opacity=opacity))
    else:
        raise NotImplementedError("Shape {} NOT YET IMPLEMENTED".format(type(shape)))


def _draw_convex_hull(svgdraw, points, convex_hull_colour, opacity):
    # plot convey hull

    last = None
    for p in _np.append(points, [points[0]], axis=0):
        if last is not None:
            l = svgdraw.line(start=last, end=p).stroke(
                width=1, color=convex_hull_colour, opacity=opacity)
            svgdraw.add(l)
        last = p

