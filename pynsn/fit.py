"""Fitting module"""
from typing import Tuple, Union
from numpy.typing import NDArray
import numpy as np
from ._stimulus import NSNStimulus
from ._shapes import PolygonShape, Dot, Ellipse, Rectangle, Picture


def total_surface_area(stim: NSNStimulus, value: float) -> None:
    """Set surface area.

    Resize all object to fit a specific surface area

    Args:
        value: surface area
    """
    stim.scale(value / stim.properties.total_surface_area)


def total_perimeter(stim: NSNStimulus, value: float) -> None:
    """fits the total parameter of the stimulus
    """
    stim.scale(value / stim.properties.total_perimeter)


def average_perimeter(stim: NSNStimulus, value: float) -> None:

    total_perimeter(stim, value * stim.n_objects)


def average_surface_area(stim: NSNStimulus, value: float) -> None:
    total_surface_area(stim, stim.n_objects * value)


# helper
def _check_polygon(stim: NSNStimulus) -> None:
    """raises exception if stimulus contains PolygonShapes"""
    if np.any(stim.shape_types == PolygonShape.name()):
        raise RuntimeError(
            "Fitting stimulus features is not possible for stimuli containing PolygonShapes")
