"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.13'

from ._lib.simple_dot_array import SimpleDotArray

from ._lib.dot_array import DotArray
from ._lib.dot_array_sequence import DASequence

from ._lib.generator import DotArrayGenerator, DASequenceGenerator, DASequenceGeneratorProcess
from ._lib.log_file import GeneratorLogger, LogFileReader

from ._lib.dot import Dot
from ._lib.colour import Colour
from ._lib.continuous_property import SurfaceArea, DotDiameter, Circumference, ConvexHull, Density
from ._lib.simple_dot_array import DotArrayProperties
from ._lib.item_features import ItemFeaturesList, ItemFeatures


# TODO:
#
#  target area
