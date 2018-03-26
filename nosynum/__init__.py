"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.6.0'

from .dot import Dot
from .dot_array import DotArrayDefinition, DotArray, list_all_saved_incremental_arrays
from .dot_array_sequences import DASequence, is_method, ALL_METHODS, M_CONVEX_HULL, M_TOTAL_AREA,\
                    M_DENSITY, M_ITEM_SIZE, M_NO_FITTING, M_TOTAL_CIRCUMFERENCE
from .multi_processing import MakeDASequenceProcess, ProcessContainer

