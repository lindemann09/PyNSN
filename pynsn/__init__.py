"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'
__version__ = '0.8'

from sys import version_info as _python_version_info
if not(_python_version_info[0] >= 3 and _python_version_info[1] >= 5):
    raise RuntimeError("PyNsN {0} ".format(__version__) +
                      "is not compatible with Python {0}.{1}. ".format(
                                                    _python_version_info[0],
                                                    _python_version_info[1]) +
                      "Please use Python 3.5+.")

from .lib._dot_array import DotArray
from .lib._colour import Colour, ImageColours
from .lib._item_attributes import ItemAttributes
from .lib._shape import Dot, Rectangle
from .lib._visual_features import Features
from .lib._match import FeatureMatcher # just for setting ITERATIVE_CONVEX_HULL_MODIFICATION and TAKE_RANDOM_DOT_FROM_CONVEXHULL toDO maybe another solution possible

from .lib import random_dot_array
from .lib import pil_image

#Further modules
# gui, pygame_surface, expyriment_stimulus,  dot_array_sequence, dot_array_archive

