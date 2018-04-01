"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

from __future__ import absolute_import

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.2'

from .lib.dot import Dot
from .lib.dot_array import DotList, DotArray, DASequence
from .lib.generator import DotArrayGenerator, DASequenceGenerator, \
                            DASequenceGeneratorProcess, GeneratorLogFile

