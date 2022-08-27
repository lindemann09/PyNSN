from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from typing import Any, Union

from numpy.typing import ArrayLike

from .._lib.coordinate import Coordinate
from ..image._colour import Colour
from .picture_file import PictureFile


class ABCShape(Coordinate, metaclass=ABCMeta):
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
        elif isinstance(attr, PictureFile):
            self._attribute = attr.attribute
        else:
            self._attribute = attr

    def get_attribute_object(self) -> Union[Colour, PictureFile, None]:
        """Class instance of the attribute, if possible

        Returns
        -------
        rtn : attribute
            If attribute represents Colour or PictureFile it returns the instance
            of the respective class otherwise None
        """

        if isinstance(self._attribute, str):
            # check if colour or picture
            col = Colour(self._attribute)
            if col.colour is not None:
                return col
            else:
                if PictureFile.is_picture_attribute(self._attribute):
                    return PictureFile(self._attribute)

        return None

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
    def rectangles_inside(self, xy: NDArray, sizes: NDArray) -> Union[np.bool_, NDArray[np.bool_]]:
        """TODO boolean or NDArray[bool] indicating whether dots are fully inside
        the shape """

    @abstractmethod
    def dots_inside(self, xy: NDArray, diameters: NDArray) -> Union[np.bool_, NDArray[np.bool_]]:
        """TODO boolean or NDArray[bool] indicating whether dots are fully inside
        the shape """
