"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.9'

from ._lib.dot import Dot
from ._lib.dot_list import DotList, DotListProperties
from ._lib.dot_array import DotArray
from ._lib.dot_array_sequence import DASequence
from ._lib.generator import DotArrayGenerator, DASequenceGenerator, DASequenceGeneratorProcess
from ._lib.files import GeneratorLogger, LogFileReader
from ._lib import colours
