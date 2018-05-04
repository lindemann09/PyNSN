"""
NonSyNum package

Creating Non-Symbolic Number Displays
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'
__version__ = '0.7.22'

from ._lib.simple_dot_array import SimpleDotArray
from ._lib.dot_array import DotArray
from ._lib.dot_array_sequence import DASequence

from ._lib.generator import DotArrayGenerator, DASequenceGenerator
from ._lib.log_file import GeneratorLogger, LogFileReader

from ._lib.dot import Dot
from ._lib.colour import Colour
from ._lib.features import TotalSurfaceArea, ItemDiameter, TotalPerimeter, FieldArea, Coverage, LogSpacing, LogSize, \
    Sparsity
from ._lib.cardinal_features import CardinalFeatures, CardinalFeaturesDotArray
from ._lib.item_attributes import ItemAttributeList, ItemAttributes

# TODO:
#
#  target area
