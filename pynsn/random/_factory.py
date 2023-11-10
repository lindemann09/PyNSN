"""
"""
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import copy
from typing import Optional, Sequence, Union, Any

from ._distributions import Levels, UnivariateDistributionType, _Constant
from .._stimulus import NSNStimulus, ShapeType
from ._appearance import Appearance


def add_random_dots(nsn_stimulus: NSNStimulus,
                    appearance: Appearance, n: int) -> NSNStimulus:
    pass


def add_random_rectangles(nsn_stimulus: NSNStimulus,
                          appearance: Appearance, n: int) -> NSNStimulus:
    pass
