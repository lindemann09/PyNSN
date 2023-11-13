"""
"""
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import copy
from typing import Optional, Sequence, Union, Any

from ._distributions import Categorical, AbstractContinuousDistr, Constant
from .._stimulus import NSNStimulus, ShapeType
from ._random_shape import AbstractRndShape


def add_random_dots(nsn_stimulus: NSNStimulus,
                    appearance: AbstractRndShape, n: int) -> NSNStimulus:
    pass


def add_random_rectangles(nsn_stimulus: NSNStimulus,
                          appearance: AbstractRndShape, n: int) -> NSNStimulus:
    pass
