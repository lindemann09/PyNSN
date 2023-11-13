"""
"""
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import copy
from typing import Optional, Sequence, Union, Any

from ._distributions import Categorical, AbstractUnivarDistr, Constant
from .._stimulus import NSNStimulus, ShapeType
from ._random_shape import AbstractRandomShape


def add_random_dots(nsn_stimulus: NSNStimulus,
                    appearance: AbstractRandomShape, n: int) -> NSNStimulus:
    pass


def add_random_rectangles(nsn_stimulus: NSNStimulus,
                          appearance: AbstractRandomShape, n: int) -> NSNStimulus:
    pass
