"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.27'

from ._lib.dot_collection import DotCollection
from ._lib.dot_array import DotArray, RandomDotArray
from ._lib.dot_array_sequence import DASequence, RandomDASequence
from ._lib.dot import Dot
from ._lib.colour import Colour
from ._lib.item_attributes import ItemAttributesList, ItemAttributes
from . import features
from ._lib.logging import LogFile, LogFileReader

# TODO:
#
#  target area
