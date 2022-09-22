"""PyNSN package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'
__version__ = '0.15.2'

from sys import version_info as _python_version_info

# pylint: disable=C0209
if not (_python_version_info[0] >= 3 and _python_version_info[1] >= 8):
    raise RuntimeError("PyNSN {0} ".format(__version__) +
                       "is not compatible with Python {0}.{1}. ".format(
                           _python_version_info[0],
                           _python_version_info[1]) +
                       "Please use Python 3.8 or later.")
# pylint: disable=C0413
from ._arrays import DotArray, RectangleArray, \
    ArrayProperties, load_array
from ._shapes import Dot, Rectangle, Picture

from ._factory.factory import NSNFactory
from ._factory import distributions

from ._lib.coordinate import Coordinate
from ._lib.colour import Colour
from ._lib.rng import init_random_generator
from ._lib import exception
from . import constants


def _print_version_info():
    from ._lib.misc import is_interactive_mode
    if is_interactive_mode():
        print("PyNSN {}".format(__version__))


_print_version_info()
