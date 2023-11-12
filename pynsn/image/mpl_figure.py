__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from typing import Optional

import numpy as _np
from matplotlib import pyplot as _plt
from matplotlib.figure import Figure as _Figure

from .. import _shapes
from .._stimulus import NSNStimulus
from . import _base
from .._stimulus.stimulus_colours import StimulusColours

# FIXME broken module


def create(
    nsn_stimulus: NSNStimulus,
    colours: Optional[StimulusColours] = None,
    dpi: float = 100,
) -> _Figure:
    """Matplotlib figure

    Parameters
    ----------
    dpi
    nsn_stimulus
    colours

    Returns
    -------

    """
    return _MatplotlibDraw().create_image(
        nsn_stimulus=nsn_stimulus, colours=colours, dpi=dpi
    )


class _MatplotlibDraw(_base.AbstractArrayDraw):
    @staticmethod
    def get_image(image_size, background_colour, **kwargs) -> _Figure:
        dpi = kwargs["dpi"]
        image_size = _np.asarray(image_size)
        figure = _plt.figure(figsize=image_size / dpi, dpi=dpi)
        if background_colour is None:
            figure.patch.set_facecolor((0, 0, 0, 0))
        else:
            figure.patch.set_facecolor(background_colour)
        axes = _plt.Axes(figure, (0, 0, 1, 1))
        axes.set_aspect("equal")  # squared
        axes.set_axis_off()
        lims = _np.transpose(image_size / 2 * _np.array([[-1, -1], [1, 1]]))
        axes.set(xlim=lims[0], ylim=lims[1])
        figure.add_axes(axes)
        return figure

    @staticmethod
    def draw_shape(img, shape: _shapes.AbstractShape, opacity):
        col = shape.colour

        if isinstance(shape, _shapes.Picture):
            raise RuntimeError(
                "Pictures are not supported for matplotlib files.")

        if isinstance(shape, _shapes.Dot):
            plt_shape = _plt.Circle(
                xy=shape.xy, radius=shape.diameter / 2, color=col.value, lw=0)
        elif isinstance(shape, _shapes.Rectangle):
            plt_shape = _plt.Rectangle(
                xy=shape.left_bottom.tolist(), width=shape.width, height=shape.height,
                color=col.value, lw=0)
        else:
            raise NotImplementedError(
                "Shape {} NOT YET IMPLEMENTED".format(type(shape)))

        plt_shape.set_alpha(opacity)
        img.axes[0].add_artist(plt_shape)

    @staticmethod
    def draw_convex_hull(img, points, convex_hull_colour, opacity):
        hull = _np.append(points, [points[0]], axis=0)
        for i in range(1, hull.shape[0]):
            line = _plt.Line2D(
                xdata=hull[i - 1: i + 1, 0],
                ydata=hull[i - 1: i + 1, 1],
                linewidth=1,
                color=convex_hull_colour.colour,
            )
            line.set_alpha(opacity)
            img.axes[0].add_artist(line)
