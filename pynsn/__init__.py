"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'
__version__ = '0.9.0'

from sys import version_info as _python_version_info
if not(_python_version_info[0] >= 3 and _python_version_info[1] >= 5):
    raise RuntimeError("PyNsN {0} ".format(__version__) +
                      "is not compatible with Python {0}.{1}. ".format(
                                                    _python_version_info[0],
                                                    _python_version_info[1]) +
                      "Please use Python 3.5+.")


from .dot_array.dot_array import DotArray
from .dot_array.item_attributes import ItemAttributes
from .dot_array.shape import Dot, Rectangle
from .dot_array.visual_features import VisualFeatures
from .dot_array.match import FeatureMatcher # just for setting ITERATIVE_CONVEX_HULL_MODIFICATION and TAKE_RANDOM_DOT_FROM_CONVEXHULL toDO maybe another solution possible
from .dot_array import random_dot_array
from .lib.colour import Colour, ImageColours
from .image import pil

#Further modules
# gui, pygame_surface, expyriment_stimulus,  dot_array_sequence, dot_array_archive

