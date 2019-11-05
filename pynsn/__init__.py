"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.38'

import sys
if not(sys.version_info[0] >= 3 and sys.version_info[1] >= 5):
    raise RuntimeError("PyNsN {0} ".format(__version__) +
                      "is not compatible with Python {0}.{1}. ".format(
                                                    sys.version_info[0],
                                                    sys.version_info[1]) +
                      "Please use Python 3.5+.")

from .lib.dot_collection import DotCollection
from .lib.dot_array import DotArray, DotArrayFactory
from .lib.dot_array_sequence import DASequence, generate_da_sequence
from .lib.geometry import Dot, Rectangle
from .lib.colour import Colour, ImageColours
from .lib.item_attributes import ItemAttributesList, ItemAttributes
from .lib import visual_features
from .lib import pil_image
from .lib.logging import LogFile, LogFileReader
