"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.5.0'

from .random_beta import random_beta, shape_parameter_beta
from . import sequence_methods
from .dot import Dot
from .dot_array import DotArrayDefinition, DotArray, list_all_saved_incremental_arrays
from .dot_array_sequences import DASequence, MakeDASequenceProcess


