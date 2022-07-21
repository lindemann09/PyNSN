__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
from matplotlib import pyplot as _plt
from .. import _shapes
from . import _array_draw


def create(object_array, colours=None, dpi=100):
    """Matplotlib figure

    Parameters
    ----------
    dpi
    object_array
    colours

    Returns
    -------

    """
    return _MatplotlibDraw().create_image(object_array=object_array, colours=colours,
                                          dpi=dpi)


class _MatplotlibDraw(_array_draw.ArrayDraw):

    @staticmethod
    def get_squared_image(image_width, background_colour, **kwargs):
        dpi = kwargs["dpi"]
        figure = _plt.figure(figsize=_np.array([image_width, image_width]) / dpi,
                             dpi=dpi)
        if background_colour is None:
            figure.patch.set_facecolor((0, 0, 0, 0))
        else:
            figure.patch.set_facecolor(background_colour)
        axes = _plt.Axes(figure, [0., 0., 1, 1])
        axes.set_aspect('equal')  # squared
        axes.set_axis_off()
        r = image_width / 2
        axes.set(xlim=[-1 * r, r], ylim=[-1 * r, r])
        figure.add_axes(axes)
        return figure

    @staticmethod
    def scale_image(figure, scaling_factor):
        # not used
        pass

    @staticmethod
    def draw_shape(img, shape, opacity, scaling_factor):
        attr = shape.get_attribute_object()

        if isinstance(attr, _shapes.PictureFile):
            raise RuntimeError("Pictures are not supported for matplotlib files.")

        if isinstance(shape, _shapes.Dot):
            r = shape.diameter / 2
            plt_shape = _plt.Circle(xy=shape.xy, radius=r, color=attr.colour,
                                    lw=0)
        elif isinstance(shape, _shapes.Rectangle):
            xy = (shape.left, shape.bottom)
            plt_shape = _plt.Rectangle(xy=xy,
                                       width=shape.width,
                                       height=shape.height,
                                       color=attr.colour, lw=0)
        else:
            raise NotImplementedError("Shape {} NOT YET IMPLEMENTED".format(type(shape)))

        plt_shape.set_alpha(opacity)
        img.axes[0].add_artist(plt_shape)

    @staticmethod
    def draw_convex_hull(img, points, convex_hull_colour, opacity,
                         scaling_factor):
        hull = _np.append(points, [points[0]], axis=0)
        for i in range(1, hull.shape[0]):
            line = _plt.Line2D(xdata=hull[i - 1:i + 1, 0],
                               ydata=hull[i - 1:i + 1, 1],
                               linewidth=1, color=convex_hull_colour.colour)
            line.set_alpha(opacity)
            img.axes[0].add_artist(line)
