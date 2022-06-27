__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
from matplotlib import pyplot as _plt
from . import _colour
from .._lib import arrays as _arrays
from .._lib import shape as _shape

def create(object_array, colours, dpi=100):
    assert isinstance(object_array, (_arrays.DotArray, _arrays.RectangleArray))
    if not isinstance(colours, _colour.ImageColours):
        raise TypeError("Colours must be of type pynsn.ImageColours")

    r = _np.ceil(object_array.target_array_radius)

    figure = _plt.figure(figsize=_np.array([r, r]) * 2 / dpi,
                         dpi=dpi)
    if colours.background.colour is None:
        figure.patch.set_facecolor((0,0,0, 0))
    else:
        figure.patch.set_facecolor(colours.background.colour)
    axes = _plt.Axes(figure, [0., 0., 1, 1])
    axes.set_aspect('equal') # squared
    axes.set_axis_off()
    axes.set(xlim=[-1*r, r], ylim=[-1*r, r])
    figure.add_axes(axes)

    if colours.target_area.colour is not None:
        obj = _shape.Dot(xy=(0, 0), diameter=object_array.target_array_radius*2,
                  attribute=colours.target_area.colour)
        _draw_shape(axes, obj)

    if object_array.features.numerosity > 0:
        if isinstance(object_array, _arrays.DotArray):
            # draw dots
            for xy, d, att in zip(object_array.xy,
                                  object_array.diameters,
                                  object_array.attributes):
                obj = _shape.Dot(xy=xy, diameter=d)
                obj.attribute = _colour.make_colour(att,
                                                    colours.default_object_colour)
                _draw_shape(axes, obj, opacity=colours.object_opacity)

        elif isinstance(object_array, _arrays.RectangleArray):
            # draw rectangle
            for xy, size, att in zip(object_array.xy,
                                     object_array.sizes,
                                     object_array.attributes):
                obj = _shape.Rectangle(xy=xy, size=size)
                obj.attribute = _colour.make_colour(att,
                                                    colours.default_object_colour)
                _draw_shape(axes, obj, opacity=colours.object_opacity)

    # draw convex hulls
    if colours.field_area_position.colour is not None:
        _draw_convex_hull(axes=axes,
                          points= object_array.features.convex_hull.position_xy,
                          convex_hull_colour=colours.field_area_position.colour,
                          opacity=colours.info_shapes_opacity)
    if colours.field_area_outer.colour is not None:
        _draw_convex_hull(axes=axes,
                          points=object_array.features.convex_hull.outer_xy,
                          convex_hull_colour=colours.field_area_outer.colour,
                          opacity=colours.info_shapes_opacity)
    #  and center of mass
    if colours.center_of_mass.colour is not None:
        obj = _shape.Dot(xy=object_array.center_of_mass(),
                         diameter=10,
                         attribute=colours.center_of_mass.colour)
        _draw_shape(axes, obj, opacity=colours.info_shapes_opacity)
    if colours.center_of_outer_positions.colour is not None:
        obj = _shape.Dot(xy=object_array.center_of_outer_positions,
                         diameter=10,
                         attribute=colours.center_of_outer_positions.colour)
        _draw_shape(axes, obj, opacity=colours.info_shapes_opacity)

    return figure


def _draw_shape(axes, shape, opacity=1.0):
    assert isinstance(shape, (_shape.Dot, _shape.Rectangle))

    colour = _colour.Colour(shape.attribute)
    if isinstance(shape, _shape.Dot):
        r = shape.diameter / 2
        plt_shape = _plt.Circle(xy=shape.xy, radius=r, color=colour.colour,
                                lw=0)
    elif isinstance(shape, _shape.Rectangle):
        xy = (shape.left, shape.bottom)
        plt_shape = _plt.Rectangle(xy=xy,
                                   width=shape.width,
                                   height=shape.height,
                                   color=colour.colour, lw=0)

    else:
        raise NotImplementedError("Shape {} NOT YET IMPLEMENTED".format(type(shape)))

    plt_shape.set_alpha(opacity)
    axes.add_artist(plt_shape)

def _draw_convex_hull(axes, points, convex_hull_colour, opacity):
    # plot convey hull
    hull = _np.append(points, [points[0]], axis=0)
    for i in range(1, hull.shape[0]):
        line = _plt.Line2D(xdata=hull[i-1:i+1, 0],
                           ydata=hull[i-1:i+1, 1],
                           linewidth = 1, color = convex_hull_colour)
        line.set_alpha(opacity)
        axes.add_artist(line)

