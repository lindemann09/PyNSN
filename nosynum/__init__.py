"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

from __future__ import absolute_import

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.3'

from .lib.dot import Dot
from .lib.dot_array import DotList, DotArray
from .lib.dot_array_sequence import DASequence
from .lib.generator import DotArrayGenerator, DASequenceGenerator, DASequenceGeneratorProcess
from .lib.files import GeneratorLogger, LogFileReader