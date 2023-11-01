"""
"""
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import copy
from typing import Optional, Sequence, Union, Any

from ._distributions import Levels, DistributionType, _Constant

ParameterDistributionType = Union[DistributionType,
                                  float, int, Sequence[Any], None]


class Appearance(object):

    def __init__(self,
                 attributes: Optional[DistributionType] = None,
                 dot_diameter: Optional[DistributionType] = None,
                 rect_width: Optional[DistributionType] = None,
                 rect_height: Optional[DistributionType] = None,
                 rect_proportion: Optional[DistributionType] = None):
        """

        Parameters
        ----------
        """

        self._attributes = _make_distr(attributes)
        self._diameter = _make_distr(dot_diameter)
        self._width = None
        self._height = None
        self._rect_proportion = None
        self.set_rectangles_sizes(width=rect_width, height=rect_height,
                                  proportion=rect_proportion)

    @property
    def attributes(self) -> Optional[DistributionType]:
        """Distribution of attributes for dots"""
        return self._attributes

    @attributes.setter
    def attributes(self, val: ParameterDistributionType):
        """Distribution of attributes for dots"""
        self._attributes = _make_distr(val)

    @property
    def dot_diameter(self) -> Optional[DistributionType]:
        """Distribution of diameter parameter"""
        return self._diameter

    @dot_diameter.setter
    def dot_diameter(self, val: ParameterDistributionType):
        """Distribution of diameter parameter"""
        self._diameter = _make_distr(val)

    @property
    def rect_width(self) -> Optional[DistributionType]:
        """Distribution of width parameter"""
        return self._width

    @property
    def rect_height(self) -> Optional[DistributionType]:
        """Distribution of height parameter"""
        return self._height

    @property
    def rect_proportion(self) -> Optional[DistributionType]:
        """Distribution of proportion parameter"""
        return self._rect_proportion

    def set_rectangles_sizes(self,
                             width: ParameterDistributionType = None,
                             height: ParameterDistributionType = None,
                             proportion: ParameterDistributionType = None):
        """Set distributions of the parameter for rectangle sizes

        Args:
            width: distribution of width (Optional)
            height: distribution of height (Optional)
            proportion: distribution of proportions (Optional)

        Notes:
            Define either rectangle width and height or rectangle proportion together
            with either width or height.

        Raises:
            TypeError: if not two of the three rectangle parameter are defined
        """
        n_rect_parameter = sum([width is not None, height is not None,
                                proportion is not None])
        if n_rect_parameter == 0:
            return
        elif n_rect_parameter != 2:
            raise TypeError("Define rectangle width and height or, alternatively, rectangle proportion together with "
                            "either width or height.")
        self._width = _make_distr(width)
        self._height = _make_distr(height)
        self._rect_proportion = _make_distr(proportion)


def _make_distr(value: ParameterDistributionType) -> Union[DistributionType, None]:
    """helper
    returns a distribution or None, if None
    """

    if value is None:
        return None
    elif isinstance(value, DistributionType):
        return value
    elif isinstance(value, (list, tuple)):
        return Levels(levels=copy(value))
    elif isinstance(value, (float, int)):
        return _Constant(value)
    else:
        raise RuntimeError("Can't make distribution from"
                           f" {type(value).__name__}")
