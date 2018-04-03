"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.5'

from ._lib.dot import Dot
from ._lib.dot_array import DotList, DotArray
from ._lib.dot_array_sequence import DASequence
from ._lib.generator import DotArrayGenerator, DASequenceGenerator, DASequenceGeneratorProcess
from ._lib.files import GeneratorLogger, LogFileReader
from ._lib import colours
