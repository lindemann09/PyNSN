"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.5.0'

from .random_beta import random_beta, shape_parameter_beta
from .dot import Dot
from .dot_array import DotArray, list_all_saved_incremental_arrays
from .dot_array_sequences import make_multiple_dot_array_sequences, make_dot_array_sequence, \
                DotArraySequence, ALL_METHODS
