from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from typing import Any, Union

from numpy.typing import ArrayLike, NDArray
from numpy import bool_

from .._lib.coordinate import Coordinate
from ..image._image_colours import Colour


class ABCShape(Coordinate, metaclass=ABCMeta):
    """Handles polar and cartesian representation (optimised processing, i.e.,
    conversions between coordinates systems will be done only once if needed)
    """

    __slots__ = ("_attribute",)

    def __init__(self, xy: ArrayLike,
                 attribute: Any) -> None:
        Coordinate.__init__(self, xy=xy)
        self._attribute = None
        self.attribute = attribute  # call setter

    @property
    def attribute(self) -> Any:
        return self._attribute

    @attribute.setter
    def attribute(self, attr: Any) -> None:
        """set attribute

        Parameters
        ----------
        attr : anything
            can be, in principle, anything.
            If Colour or PictureFile, it will convert it to their string
            representation
        """
        if isinstance(attr, Colour):
            self._attribute = attr.colour
        else:
            self._attribute = attr

    def get_colour(self) -> Colour:
        """Class instance of the attribute, if possible

        Returns
        -------
        rtn : Colour
        """

        if isinstance(self._attribute, str):
            # check if colour or picture
            col = Colour(self._attribute)
            if col.colour is not None:
                return col
        return Colour(None)

    @abstractmethod
    def __repr__(self) -> str:
        """"""

    @abstractmethod
    def distance(self, other: ABCShape) -> float:
        """Euclidean distance to another shapes

        Notes:
            Negative distances indicate an overlap and represent minimum distance
            between outer shape boarders
        """

    @property
    @abstractmethod
    def area(self) -> float:
        """area of the object/shape"""

    @property
    @abstractmethod
    def perimeter(self) -> float:
        """perimeter"""

    @abstractmethod
    def rectangles_inside(self, xy: NDArray, sizes: NDArray) -> Union[bool_, NDArray[bool_]]:
        """TODO boolean or NDArray[bool] indicating whether dots are fully inside
        the shape """

    @abstractmethod
    def dots_inside(self, xy: NDArray, diameters: NDArray) -> Union[bool_, NDArray[bool_]]:
        """TODO boolean or NDArray[bool] indicating whether dots are fully inside
        the shape """
