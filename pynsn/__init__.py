"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'
__version__ = '0.11.6'

from sys import version_info as _python_version_info
if not(_python_version_info[0] >= 3 and _python_version_info[1] >= 6):
    raise RuntimeError("PyNsN {0} ".format(__version__) +
                      "is not compatible with Python {0}.{1}. ".format(
                                                    _python_version_info[0],
                                                    _python_version_info[1]) +
                      "Please use Python 3.6 or later.")

from ._lib import Point, Dot, Rectangle, \
                AttributeArray, DotArray, RectangleArray, \
                NoSolutionError, PictureFile

from .image._colour import Colour, ImageColours

from . import random_array
from . import visual_properties
from . import distributions

