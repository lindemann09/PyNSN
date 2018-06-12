"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.26'

from ._lib.simple_dot_array import SimpleDotArray
from ._lib.dot_array import DotArray
from ._lib.dot_array_sequence import DASequence
from ._lib.generator import DotArrayGenerator, DASequenceGenerator
from ._lib.dot import Dot
from ._lib.colour import Colour
from ._lib.item_attributes import ItemAttributeList, ItemAttributes
from . import features
from ._lib.logging import LogFile, LogFileReader

# TODO:
#
#  target area
