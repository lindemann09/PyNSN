__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as _np
from matplotlib import pyplot as _pyplot
from . import _colour
from .._lib  import arrays as _arrays
from .._lib.geometry import cartesian2image_coordinates as _c2i_pos

def create(object_array, colours):
    assert isinstance(object_array, (_arrays.DotArray, _arrays.RectangleArray))
    assert isinstance(colours, _colour.ImageColours)

    raise NotImplementedError()