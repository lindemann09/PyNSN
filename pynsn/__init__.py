"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.36'

import sys
if not(sys.version_info[0] >= 3 and sys.version_info[1] >= 5):
    raise RuntimeError("PyNsN {0} ".format(__version__) +
                      "is not compatible with Python {0}.{1}. ".format(
                                                    sys.version_info[0],
                                                    sys.version_info[1]) +
                      "Please use Python 3.5+.")

from ._lib.dot_collection import DotCollection
from ._lib.dot_array import DotArray, DotArrayGenerator
from ._lib.dot_array_sequence import DASequence, generate_da_sequence
from ._lib.geometry import Dot, Rectangle
from ._lib.colour import Colour
from ._lib.item_attributes import ItemAttributesList, ItemAttributes
from ._lib import visual_features
from ._lib.logging import LogFile, LogFileReader
from ._lib.pil_image import PILImagePlotter