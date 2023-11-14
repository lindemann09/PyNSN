"""Fitting module"""
from typing import Tuple, Union
from numpy.typing import NDArray
import numpy as np
from ._stimulus import NSNStimulus
from ._shapes import PolygonShape, Dot, Ellipse, Rectangle, Picture



def scale_size(stim:NSNStimulus, scaling:Union[float, np.float_]) -> None:
    """scales the size of the shapes in the stimulus"""
    _check_polygon(stim) #only rectangles and circular shapes only

    # pylint: disable=W0212
    stim._sizes= stim._sizes* scaling
    stim._reset()


def total_surface_area(stim:NSNStimulus, value: float) -> None:
    """Set surface area.

    Resize all object to fit a specific surface area

    Args:
        value: surface area
    """

    scale_size(stim, scaling=value / stim.properties.total_surface_area)


def total_perimeter(stim:NSNStimulus, value: float) -> None:
    """fits the total parameter of the stimulus
    """
    scale_size(stim, scaling=value / stim.properties.total_perimeter)

def average_perimeter(stim:NSNStimulus, value: float) -> None:

    total_perimeter(stim, value * stim.n_objects)

def average_surface_area(stim:NSNStimulus, value: float) -> None:
    total_surface_area(stim, stim.n_objects * value)



# helper
def _check_polygon(stim:NSNStimulus) -> None:
    """raises exception if stimulus contains PolygonShapes"""
    if np.any(stim.shape_types == PolygonShape.name()):
        raise RuntimeError("Fitting stimulus features is not possible for stimuli containing PolygonShapes")
