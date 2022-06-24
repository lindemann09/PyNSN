__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
from matplotlib import pyplot as _plt
from . import _colour
from .._lib import arrays as _arrays
from .._lib import shape as _shape

def create(object_array, colours):
    assert isinstance(object_array, (_arrays.DotArray, _arrays.RectangleArray))
    if not isinstance(colours, _colour.ImageColours):
        raise TypeError("Colours must be of type pynsn.ImageColours")

    figure, axes = _plt.subplots()
    lims = _np.ceil([-object_array.target_array_radius, object_array.target_array_radius])
    axes.set_aspect(1) # squared
    axes.set(xlim=lims, ylim=lims)

    if colours.target_area.colour is not None:
        axes.add_artist(_plt.Circle(xy=(0, 0),
                                    radius=object_array.target_array_radius,
                                    color=colours.target_area.colour))
    #FIXME Size dots too big (only plot filled)
    #FIXME Rectangle not implemented

    for xy, d, att in zip(object_array.xy,
                          object_array.diameters,
                          object_array.attributes):
        if att is None:
            c = colours.default_item_colour
        else:
            try:
                c = _colour.Colour(att)
            except TypeError:
                c = colours.default_item_colour

        axes.add_artist(_plt.Circle(xy=xy, radius=d/2, color=c.colour))

    return figure, axes

