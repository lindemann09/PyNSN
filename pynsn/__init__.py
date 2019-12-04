"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.47'

from sys import version_info as _python_version_info
if not(_python_version_info[0] >= 3 and _python_version_info[1] >= 5):
    raise RuntimeError("PyNsN {0} ".format(__version__) +
                      "is not compatible with Python {0}.{1}. ".format(
                                                    _python_version_info[0],
                                                    _python_version_info[1]) +
                      "Please use Python 3.5+.")

from ._lib._dot_array import DotArray
from ._lib._colour import Colour, ImageColours
from ._lib._item_attributes import ItemAttributes
from ._lib._shape import Dot, Rectangle
from ._lib._visual_features import Features
from ._lib._match import FeatureMatcher # just for setting ITERATIVE_CONVEX_HULL_MODIFICATION and TAKE_RANDOM_DOT_FROM_CONVEXHULL toDO maybe another solution possible

from ._lib import random_dot_array
from ._lib import pil_image

#Further modules
# gui, pygame_surface, expyriment_stimulus,  dot_array_sequence, dot_array_archive

