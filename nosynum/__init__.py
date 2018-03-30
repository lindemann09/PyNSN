"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.0'

from .dot import Dot
from .dot_array import NumpyDotList, DotArray
from .dot_array_sequences import DASequence, is_method, ALL_METHODS, M_CONVEX_HULL, M_TOTAL_AREA,\
                    M_DENSITY, M_ITEM_SIZE, M_NO_FITTING, M_TOTAL_CIRCUMFERENCE
from .random_da_creator import RandomDotArrayCreater
from .multi_processing import MakeDASequenceProcess, ProcessContainer

