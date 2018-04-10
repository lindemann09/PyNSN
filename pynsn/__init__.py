"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.10'

from ._lib.dot import Dot
from ._lib.simple_dot_array import SimpleDotArray, DotArrayProperties
from ._lib.dot_array import DotArray
from ._lib.dot_array_sequence import DASequence
from ._lib.generator import DotArrayGenerator, DASequenceGenerator, DASequenceGeneratorProcess
from ._lib.files import GeneratorLogger, LogFileReader
from ._lib import colours


# TODO:
#
#  * other shapes
# * rectangluar array
