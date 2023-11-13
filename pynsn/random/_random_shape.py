"""
"""
from pathlib import Path

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from typing import Optional, Sequence, Union

from numpy.typing import NDArray

from ._distributions import (AbstractUnivarDistr, Categorical, CategoricalLike,
                             Constant, ConstantLike)

DistributionLike = Union[AbstractUnivarDistr, ConstantLike, CategoricalLike]


def _make_distr(value: DistributionLike) -> AbstractUnivarDistr:
    """helper"""

    if isinstance(value, AbstractUnivarDistr):
        return value
    elif isinstance(value, ConstantLike): # constant like
        return Constant(value)
    else:
        return Categorical(value)


class AbstractRndShape(metaclass=ABCMeta):

    def __init__(self, attributes: Optional[DistributionLike] = None):
        """

        Parameters
        ----------
        """
        if attributes is None:
            self._attributes = None
        else:
            self._attributes = _make_distr(attributes)

    @property
    def attributes(self) -> Optional[AbstractUnivarDistr]:
        """Distribution of attributes """
        return self._attributes

    @classmethod
    @property
    def name(cls) -> str:
        return str(cls.__name__)

    def __repr__(self) -> str:
        d = self.todict()
        del d['type']
        return f"{self.name}({d})"

    @abstractmethod
    def todict(self) -> dict:
        """dict representation of the object"""
        if  isinstance(self.attributes, AbstractUnivarDistr):
            attr = self.attributes.todict()
        elif self.attributes is None:
            attr = None
        else:
            attr = str(self.attributes)
        return {"type": self.name, "attr": attr}


class RndDot(AbstractRndShape):

    def __init__(self,
                 diameter: DistributionLike,
                 attributes: Optional[DistributionLike] = None):
        """Define distributions parameter
        """
        super().__init__(attributes=attributes)
        self._diameter = _make_distr(diameter)

    @property
    def diameter(self) -> Optional[AbstractUnivarDistr]:
        return self._diameter

    def todict(self) -> dict:
        rtn = super().todict()
        if isinstance(self._diameter, AbstractUnivarDistr):
            d = self._diameter.todict()
        else:
            d = None
        rtn.update({"diameter": d})
        return rtn

class _RandShapeWidthHeight(AbstractRndShape, metaclass=ABCMeta):

    def __init__(self,
                 width: Optional[DistributionLike] = None,
                 height: Optional[DistributionLike] = None,
                 size_proportion: Optional[DistributionLike] = None,
                 attributes: Optional[DistributionLike] = None):
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
            self._width = _make_distr(width)
        if height is None:
            self._height = None
        else:
            self._height = _make_distr(height)
        if size_proportion is None:
            self._size_proportion = None
        else:
            self._size_proportion = _make_distr(size_proportion)

    @property
    def width(self) -> Optional[AbstractUnivarDistr]:
        """Distribution of width parameter"""
        return self._width

    @property
    def height(self) -> Optional[AbstractUnivarDistr]:
        """Distribution of height parameter"""
        return self._height

    @property
    def size_proportion(self) -> Optional[AbstractUnivarDistr]:
        """Distribution of proportion parameter"""
        return self._size_proportion

    def todict(self) -> dict:
        rtn = super().todict()
        if self._width is None:
            w = None
        else:
            w = self._width.todict()
        if self._height is None:
            h = None
        else:
            h = self._height.todict()
        if self._size_proportion is None:
            s = None
        else:
            s = self._size_proportion.todict()
        rtn.update({"width": w, "height":h, "size_proportion": s})
        return rtn


class RndRectangle(_RandShapeWidthHeight):
    pass


class RndEllipse(_RandShapeWidthHeight):
    pass

class RndPolygon(_RandShapeWidthHeight):
    pass

class RndPicture(_RandShapeWidthHeight):

    def __init__(self,
                 width: Optional[DistributionLike] = None,
                 height: Optional[DistributionLike] = None,
                 size_proportion: Optional[DistributionLike] = None,
                 path: Union[Categorical, Constant, str, Path, Sequence[Path],
                             Sequence[str], None] = None):
        """Define distributions parameter
        """
        if isinstance(path, (str, Path)):
            attr = Constant(str(Path(path)))
        elif isinstance(path, Sequence):
            pathes = [Path(x)for x in path]
            attr = Categorical(pathes)
        else:
            attr = path
        super().__init__(width=width, height=height, size_proportion=size_proportion,
                         attributes=attr)


