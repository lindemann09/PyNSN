""" """

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from copy import deepcopy
from pathlib import Path
from typing import List, Sequence

import numpy as np
import shapely
from numpy.typing import NDArray

from .._shapes import AbstractShape, Dot, Ellipse, Picture, PolygonShape, Rectangle
from ._distributions import (
    AbstractUnivarDistr,
    Categorical,
    Constant,
    ConstantLike,
    DistributionLike,
)


def _make_distr(value: DistributionLike) -> AbstractUnivarDistr:
    """helper"""

    if isinstance(value, AbstractUnivarDistr):
        return value
    elif isinstance(value, ConstantLike):  # constant like
        return Constant(value)
    else:
        return Categorical(value)


class AbstractRndShape(metaclass=ABCMeta):
    def __init__(self, attributes: DistributionLike | None = None):
        """

        Parameters
        ----------
        """
        if attributes is None:
            self._attributes = None
        else:
            self._attributes = _make_distr(attributes)

    @property
    def attributes(self) -> AbstractUnivarDistr | None:
        """Distribution of attributes"""
        return self._attributes

    @classmethod
    def shape_type(cls) -> str:
        return cls.__name__

    def __repr__(self) -> str:
        d = self.to_dict()
        del d["type"]
        return f"{self.shape_type()}({d})"

    @abstractmethod
    def to_dict(self) -> dict:
        """dict representation of the object"""
        if isinstance(self.attributes, AbstractUnivarDistr):
            attr = self.attributes.to_dict()
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
    def __init__(
        self, diameter: DistributionLike, attributes: DistributionLike | None = None
    ):
        """Define distributions parameter"""
        super().__init__(attributes=attributes)
        self._diameter = _make_distr(diameter)

    @property
    def diameter(self) -> AbstractUnivarDistr | None:
        return self._diameter

    def to_dict(self) -> dict:
        rtn = super().to_dict()
        if isinstance(self._diameter, AbstractUnivarDistr):
            d = self._diameter.to_dict()
        else:
            d = None
        rtn.update({"diameter": d})
        return rtn

    def sample(self, n: int = 1) -> List[Dot]:
        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)
        return [
            Dot(diameter=d, attribute=a) for d, a in zip(self._diameter.sample(n), attr)
        ]


class _RandShapeWidthHeight(AbstractRndShape, metaclass=ABCMeta):
    def __init__(
        self,
        width: DistributionLike | None = None,
        height: DistributionLike | None = None,
        size_proportion: DistributionLike | None = None,
        attributes: DistributionLike | None = None,
    ):
        """Define distributions parameter

        Parameters
        ----------
        width : DistributionLike | None, optional
            distribution of width , by default None
        height : DistributionLike | None, optional
            distribution of height, by default None
        size_proportion : DistributionLike | None, optional
            distribution of proportions width/height, by default None
        attributes : DistributionLike | None, optional
            distribution of attributes, by default None

        Raises
        ------
        TypeError
            if not two of the three rectangle parameter are defined

        Notes
        -----
        Define either rectangle width and height or rectangle proportion together
        with either width or height.

        """

        super().__init__(attributes)
        n_parameter = sum([width is not None, height is not None])
        if size_proportion is not None:
            if n_parameter != 1:
                raise TypeError(
                    "Define size proportion together with either width or height, not both."
                )
        elif n_parameter < 2:
            raise TypeError(
                "Define width and height or, alternatively, size proportion together with width or height."
            )

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
    def width(self) -> AbstractUnivarDistr | None:
        """Distribution of width parameter"""
        return self._width

    @property
    def height(self) -> AbstractUnivarDistr | None:
        """Distribution of height parameter"""
        return self._height

    @property
    def size_proportion(self) -> AbstractUnivarDistr | None:
        """Distribution of proportion parameter (width/height)"""
        return self._size_proportion

    def to_dict(self) -> dict:
        rtn = super().to_dict()
        if self._width is None:
            w = None
        else:
            w = self._width.to_dict()
        if self._height is None:
            h = None
        else:
            h = self._height.to_dict()
        if self._size_proportion is None:
            s = None
        else:
            s = self._size_proportion.to_dict()
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
            raise RuntimeError(
                "Something went wrong with the definition of the" " RndRectangle"
            )
        return size


class RndRectangle(_RandShapeWidthHeight):
    def sample(self, n: int = 1) -> List[Rectangle]:
        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)
        return [
            Rectangle(size=s, attribute=a) for s, a in zip(self._sample_sizes(n), attr)
        ]


class RndEllipse(_RandShapeWidthHeight):
    def sample(self, n: int = 1) -> List[Ellipse]:
        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)
        size = self._sample_sizes(n)
        return [Ellipse(size=s, attribute=a) for s, a in zip(size, attr)]


class RndPolygonShape(_RandShapeWidthHeight):
    def __init__(
        self,
        polygons: (shapely.Polygon | Sequence[shapely.Polygon] | NDArray),  # type: ignore
        width: DistributionLike | None = None,
        height: DistributionLike | None = None,
        size_proportion: DistributionLike | None = None,
        attributes: DistributionLike | None = None,
    ):
        """Define distributions parameter


        Parameters
        ----------
        polygons : shapely.Polygon  |  Sequence[shapely.Polygon]  |  NDArray[shapely.Polygon]
            _description_
        width : DistributionLike | None, optional
            distribution of width, by default None
        height : DistributionLike | None, optional
            distribution of height, by default None
        size_proportion : DistributionLike | None, optional
            distribution of proportions width/height , by default None
        attributes : DistributionLike | None, optional
            distribution of attributes, by default None

        Raises
        ------
        TypeError :
            if not two of the three rectangle parameter are defined

        Notes
        ------
            Define either rectangle width and height or rectangle proportion together
            with either width or height.

        """
        if width is None and height is None:
            width = -1
            height = -1
            size_proportion = None
        super().__init__(
            width=width,
            height=height,
            size_proportion=size_proportion,
            attributes=attributes,
        )

        if isinstance(polygons, shapely.Polygon):
            self._polygons = Constant(PolygonShape(polygons))
        else:
            self._polygons = Categorical([PolygonShape(p) for p in polygons])

    @property
    def polygons(self) -> Constant | Categorical:
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
    def __init__(
        self,
        width: DistributionLike | None = None,
        height: DistributionLike | None = None,
        size_proportion: DistributionLike | None = None,
        path: Categorical
        | Constant
        | str
        | Path
        | Sequence[Path]
        | Sequence[str]
        | None = None,
    ):
        """Define distributions parameter"""
        if isinstance(path, (str, Path)):
            attr = Constant(str(Path(path)))
        elif isinstance(path, Sequence):
            pathes = [Path(x) for x in path]
            attr = Categorical(pathes)
        else:
            attr = path
        super().__init__(
            width=width, height=height, size_proportion=size_proportion, attributes=attr
        )

    def sample(self, n: int = 1) -> List[Picture]:
        if self._attributes is None:
            attr = [None] * n
        else:
            attr = self._attributes.sample(n)
        return [Picture(size=s, path=a) for s, a in zip(self._sample_sizes(n), attr)]  # type: ignore
