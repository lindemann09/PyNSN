"""
"""

from copy import deepcopy
import shapely
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import List, Optional, Sequence, Union

import numpy as np
from numpy.typing import NDArray

from ._distributions import (AbstractUnivarDistr, Categorical, CategoricalLike,
                             Constant, ConstantLike)

from .._shapes import Dot, Ellipse, Rectangle, Picture, PolygonShape, AbstractShape

DistributionLike = Union[AbstractUnivarDistr, ConstantLike, CategoricalLike]


def _make_distr(value: DistributionLike) -> AbstractUnivarDistr:
    """helper"""

    if isinstance(value, AbstractUnivarDistr):
        return value
    elif isinstance(value, ConstantLike):  # constant like
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
    def shape_type(cls) -> str:
        return cls.__name__

    def __repr__(self) -> str:
        d = self.todict()
        del d['type']
        return f"{self.shape_type()}({d})"

    @abstractmethod
    def todict(self) -> dict:
        """dict representation of the object"""
        if isinstance(self.attributes, AbstractUnivarDistr):
            attr = self.attributes.todict()
        elif self.attributes is None:
            attr = None
        else:
            attr = str(self.attributes)
        return {"type": self.shape_type(), "attr": attr}

    @abstractmethod
    def sample(self, n: int = 1) -> List[AbstractShape]:
        """get n random variants of shapes

        Note
        -----
        position is xy=(0,0), except for PolygonShapes the position is the
        centroid of the shapley.polygon.
        """


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

    def sample(self, n: int = 1) -> List[Dot]:

        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)
        return [Dot(diameter=d, attribute=a)
                for d, a in zip(self._diameter.sample(n), attr)]


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
            proportion: distribution of proportions width/height (Optional)

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
        """Distribution of proportion parameter (width/height)"""
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
        rtn.update({"width": w, "height": h, "size_proportion": s})
        return rtn

    def _sample_sizes(self, n: int) -> NDArray:

        size = np.empty((n, 2), dtype=float)
        if self._width is not None and self._height is not None:
            size[:, 0] = self._width.sample(n)
            size[:, 1] = self._height.sample(n)
        elif self._width is not None and self._size_proportion is not None:
            size[:, 0] = self._width.sample(n)
            size[:, 1] = size[:, 0] / self._size_proportion.sample(n)
        elif self._height is not None and self._size_proportion is not None:
            size[:, 1] = self._height.sample(n)
            size[:, 0] = size[:, 1] * self._size_proportion.sample(n)
        else:
            raise RuntimeError("Something went wrong with the definition of the"
                               " RndRectangle")
        return size


class RndRectangle(_RandShapeWidthHeight):

    def sample(self, n: int = 1) -> List[Rectangle]:
        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)
        return [Rectangle(size=s, attribute=a)
                for s, a in zip(self._sample_sizes(n), attr)]


class RndEllipse(_RandShapeWidthHeight):

    def sample(self, n: int = 1) -> List[Ellipse]:
        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)
        size = self._sample_sizes(n)
        return [Ellipse(size=s, attribute=a) for s, a in zip(size, attr)]


class RndPolygonShape(_RandShapeWidthHeight):

    def __init__(self,
                 # type: ignore
                 polygons: Union[shapely.Polygon, Sequence[shapely.Polygon], NDArray[shapely.Polygon]],
                 width: Optional[DistributionLike] = None,
                 height: Optional[DistributionLike] = None,
                 size_proportion: Optional[DistributionLike] = None,
                 attributes: Optional[DistributionLike] = None):
        """Define distributions parameter

        Args:
            width: distribution of width (Optional)
            height: distribution of height (Optional)
            proportion: distribution of proportions width/height (Optional)

        Notes:
            Define either rectangle width and height or rectangle proportion together
            with either width or height.

        Raises:
            TypeError: if not two of the three rectangle parameter are defined
        """
        if width is None and height is None:
            width = -1
            height = -1
            size_proportion = None
        super().__init__(width=width, height=height,
                         size_proportion=size_proportion, attributes=attributes)

        if isinstance(polygons, shapely.Polygon):
            self._polygons = Constant(PolygonShape(polygons))
        else:
            self._polygons = Categorical([PolygonShape(p) for p in polygons])

    @property
    def polygons(self) -> Union[Constant, Categorical]:
        return self._polygons

    def sample(self, n: int = 1) -> List[PolygonShape]:
        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)

        rtn = []
        for s, p, a in zip(self._sample_sizes(n), self._polygons.sample(n), attr):
            assert isinstance(p, PolygonShape)
            if s[0] > 0:
                x = deepcopy(p)
                x.size = s
                x.attribute = a
                rtn.append(x)
        return rtn


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

    def sample(self, n: int = 1) -> List[Picture]:
        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)
        return [Picture(size=s,
                        path=a) for s, a in zip(self._sample_sizes(n), attr)]  # type: ignore
