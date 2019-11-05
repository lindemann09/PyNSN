"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.39'

from sys import version_info as _version_info
if not(_version_info[0] >= 3 and _version_info[1] >= 5):
    raise RuntimeError("PyNsN {0} ".format(__version__) +
                      "is not compatible with Python {0}.{1}. ".format(
                                                    _version_info[0],
                                                    _version_info[1]) +
                      "Please use Python 3.5+.")

from ._lib._dot_array import DotCollection, DotArray
from ._lib._colour import Colour, ImageColours
from ._lib._item_attributes import ItemAttributesList, ItemAttributes
from ._lib._logging import LogFile, LogFileReader

from ._lib import visual_features
from ._lib import pil_image
from ._lib import geometry
from ._lib import random_dot_array
from ._lib import dot_array_sequence
