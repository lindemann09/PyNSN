__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
from matplotlib import pyplot as _plt
from . import _colour
from .. import _lib

# FIXME can't handle pictures


def create(object_array, colours=None, dpi=100):
    """create a matplotlib figure"""

    _lib._check_object_array(object_array)
    if colours is None:
        colours = _colour.ImageColours()
    if not isinstance(colours, _colour.ImageColours):
        raise TypeError("Colours must be of type image.ImageColours")

    r = _np.ceil(object_array.target_area_radius)

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
        obj = _lib.Dot(xy=(0, 0), diameter=object_array.target_area_radius * 2,
                          attribute=colours.target_area.colour)
        _draw_shape(axes, obj)

    if object_array.properties.numerosity > 0:
        if isinstance(object_array, _lib.DotArray):
            # draw dots
            for xy, d, att in zip(object_array.xy,
                                  object_array.diameters,
                                  object_array.attributes):
                obj = _lib.Dot(xy=xy, diameter=d)
                obj.attribute = _colour.Colour(att,
                                                    colours.default_object_colour)
                _draw_shape(axes, obj, opacity=colours.opacity_object)

        elif isinstance(object_array, _lib.RectangleArray):
            # draw rectangle
            for xy, size, att in zip(object_array.xy,
                                     object_array.sizes,
                                     object_array.attributes):
                obj = _lib.Rectangle(xy=xy, size=size)
                obj.attribute = _colour.Colour(att,
                                        colours.default_object_colour)
                _draw_shape(axes, obj, opacity=colours.opacity_object)

    # draw convex hulls
    if colours.field_area_positions.colour is not None and \
            object_array.properties.field_area_positions > 0:
        _draw_convex_hull(axes=axes,
                          points= object_array.properties.convex_hull_positions.xy,
                          convex_hull_colour=colours.field_area_positions.colour,
                          opacity=colours.opacity_guides)
    if colours.field_area.colour is not None and \
            object_array.properties.field_area > 0:
        _draw_convex_hull(axes=axes,
                          points=object_array.properties.convex_hull.xy,
                          convex_hull_colour=colours.field_area.colour,
                          opacity=colours.opacity_guides)
    #  and center of mass
    if colours.center_of_field_area.colour is not None:
        obj = _lib.Dot(xy=object_array.center_of_field_area(),
                          diameter=10,
                          attribute=colours.center_of_field_area.colour)
        _draw_shape(axes, obj, opacity=colours.opacity_guides)
    if colours.center_of_mass.colour is not None:
        obj = _lib.Dot(xy=object_array.center_of_mass(),
                          diameter=10,
                          attribute=colours.center_of_mass.colour)
        _draw_shape(axes, obj, opacity=colours.opacity_guides)

    return figure


def _draw_shape(axes, shape, opacity=1.0):
    assert isinstance(shape, (_lib.Dot, _lib.Rectangle))

    colour = _colour.Colour(shape.attribute)
    if isinstance(shape, _lib.Dot):
        r = shape.diameter / 2
        plt_shape = _plt.Circle(xy=shape.xy, radius=r, color=colour.colour,
                                lw=0)
    elif isinstance(shape, _lib.Rectangle):
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

