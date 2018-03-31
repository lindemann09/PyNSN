"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

from __future__ import absolute_import

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.0'

from .lib.dot import Dot
from .lib.dot_array import DotArray, NumpyDotList
from .lib.dot_array_sequences import DASequence, is_method, ALL_METHODS, M_CONVEX_HULL, M_TOTAL_AREA,\
                    M_DENSITY, M_ITEM_SIZE, M_NO_FITTING, M_TOTAL_CIRCUMFERENCE
from .lib.generator import RandomDotArrayGenerator
from .lib.multi_processing import MakeDASequenceProcess, ProcessContainer

