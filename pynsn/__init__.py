"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'
__version__ = '0.9.5-3'

from sys import version_info as _python_version_info
if not(_python_version_info[0] >= 3 and _python_version_info[1] >= 5):
    raise RuntimeError("PyNsN {0} ".format(__version__) +
                      "is not compatible with Python {0}.{1}. ".format(
                                                    _python_version_info[0],
                                                    _python_version_info[1]) +
                      "Please use Python 3.5+.")


from .nsn.shape import Dot, Rectangle
from .nsn.colour import Colour, ImageColours
from .nsn.dot_array import DotArray
from .nsn import factory

from .nsn.visual_features import VisualFeatures
#Further modules
# gui, pygame_surface, expyriment_stimulus,  dot_array_sequence, dot_array_archive

