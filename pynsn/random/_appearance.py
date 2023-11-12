"""
"""
from pathlib import Path
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from numpy.typing import NDArray
from abc import ABCMeta
from typing import Optional, Sequence, Tuple, Union, Any
from .._shapes.colour import Colour
from ._distributions import Categorical, AbstractUnivarDistr, Constant, AbstractDistribution

ConstantLike = Union[Constant, float, int, str, Colour]
UnivariateLike = Union[ConstantLike, AbstractUnivarDistr]
AttributesDistrLike = Union[ConstantLike,
                            AbstractDistribution, Sequence, NDArray, None]


def _make_univar_distr(value: UnivariateLike) -> Union[Constant, AbstractUnivarDistr]:
    """helper"""

    if isinstance(value, AbstractUnivarDistr):
        return value
    else:  # constant like
        return Constant(value)


class Appearance(metaclass=ABCMeta):

    def __init__(self, attributes: AttributesDistrLike = None):
        """

        Parameters
        ----------
        """
        if attributes is None or isinstance(attributes, AbstractDistribution):
            self._attributes = None
        elif isinstance(attributes, ConstantLike):
            self._attributes = Constant(attributes)
        else:
            self._attributes = Categorical(categories=attributes)

    @property
    def attributes(self) -> Optional[AbstractDistribution]:
        """Distribution of attributes """
        return self._attributes


class RandomDot(Appearance):

    def __init__(self,
                 diameter: UnivariateLike,
                 attributes: AttributesDistrLike = None):
        """Define distributions parameter
        """
        super().__init__(attributes=attributes)
        self._diameter = _make_univar_distr(diameter)

    @property
    def diameter(self) -> Union[Constant, AbstractUnivarDistr]:
        return self._diameter


class _AppearanceWidthHeight(Appearance, metaclass=ABCMeta):

    def __init__(self,
                 width: Optional[UnivariateLike] = None,
                 height: Optional[UnivariateLike] = None,
                 size_proportion: Optional[UnivariateLike] = None,
                 attributes: AttributesDistrLike = None):
        """Define distributions parameter

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

        super().__init__(attributes)
        n_parameter = sum([width is not None, height is not None])
        if size_proportion is not None:
            if n_parameter != 1:
                raise TypeError(
                    "Define size proportion together with either width or height, not both.")
        elif n_parameter < 2:
            raise TypeError(
                "Define width and height or, alternatively, size proportion together with width or height.")

        if width is None:
            self._width = None
        else:
            self._width = _make_univar_distr(width)
        if height is None:
            self._height = None
        else:
            self._height = _make_univar_distr(height)
        if size_proportion is None:
            self._size_proportion = None
        else:
            self._size_proportion = _make_univar_distr(size_proportion)

    @property
    def width(self) -> Union[Constant, AbstractUnivarDistr, None]:
        """Distribution of width parameter"""
        return self._width

    @property
    def height(self) -> Union[Constant, AbstractUnivarDistr, None]:
        """Distribution of height parameter"""
        return self._height

    @property
    def size_proportion(self) -> Union[Constant, AbstractUnivarDistr, None]:
        """Distribution of proportion parameter"""
        return self._size_proportion


class RandomRectangle(_AppearanceWidthHeight):
    pass


class RandomEllipse(_AppearanceWidthHeight):
    pass


class RandomPicture(_AppearanceWidthHeight):

    def __init__(self,
                 width: Optional[UnivariateLike] = None,
                 height: Optional[UnivariateLike] = None,
                 size_proportion: Optional[UnivariateLike] = None,
                 path: Union[Categorical, Constant, str, Path, Sequence[Path],
                             Sequence[str], None] = None):
        """Define distributions parameter
        """
        if isinstance(path, (str, Path)):
            attr = Constant(Path(path))
        elif isinstance(path, Sequence):
            pathes = [Path(x)for x in path]
            attr = Categorical(pathes)
        else:
            attr = path
        super().__init__(width=width, height=height, size_proportion=size_proportion,
                         attributes=attr)
