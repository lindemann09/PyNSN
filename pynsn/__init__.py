"""
PyNSN package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'
__version__ = '0.12.1'


from sys import version_info as _python_version_info

if not (_python_version_info[0] >= 3 and _python_version_info[1] >= 8):
    raise RuntimeError("PyNSN {0} ".format(__version__) +
                       "is not compatible with Python {0}.{1}. ".format(
                           _python_version_info[0],
                           _python_version_info[1]) +
                       "Please use Python 3.8 or later.")


from ._lib import Coordinate, Point, Dot, Rectangle, \
    DotArray, RectangleArray, PointArray, ArrayParameter, \
    PictureFile, NSNFactory

from .image._colour import Colour, ImageColours
from ._lib.rng import init_random_generator

from . import visual_properties
from . import distributions


def _print_version_info():
    from ._lib.misc import is_interactive_mode
    if is_interactive_mode():
        print("PyNSN {}".format(__version__))

_print_version_info()
