"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

from __future__ import absolute_import

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.1'

from .lib.dot import Dot
from .lib.dot_array import DotArray, NumpyDotList, DotArraySequence
from .lib.generator import RandomDotArrayGenerator, DASequenceGenerator, DASequenceGeneratorProcess,\
                    GeneratorLogFile,                    \
                    ALL_METHODS, M_CONVEX_HULL, M_TOTAL_AREA,\
                    M_DENSITY, M_MEAN_DIAMETER, M_NO_FITTING, M_TOTAL_CIRCUMFERENCE

