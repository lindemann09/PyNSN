from abc import ABCMeta, abstractmethod

from typing import Any, Dict, List, Optional, Iterator
from numpy.typing import ArrayLike, NDArray
from .._lib.typing import IntOVector


class ABCObjectArray(metaclass=ABCMeta):
    """Class of numpy representation of dot"""

    @abstractmethod
    def check(self) -> None:
        """raises value error if badly shaped data"""

    @abstractmethod
    def set_attributes(self, attributes: Optional[ArrayLike]) -> None:
        """Set all attributes

        Args:
            attributes: single attribute or list of attributes
        """

    @abstractmethod
    def np_append(self,
                  xy: ArrayLike,
                  diameter: ArrayLike,
                  attributes: Optional[ArrayLike] = None) -> None:
        """if attributes comprises only one element all new objects will get
        this attribute """

    @abstractmethod
    def delete(self, index: IntOVector) -> None:
        """ """

    @abstractmethod
    def clear(self):
        """ """

    @abstractmethod
    def add(self, shapes) -> None:
        """append one dot or list of dots"""

    @property
    @abstractmethod
    def surface_areas(self) -> NDArray:
        """surface area of each dot"""

    @property
    @abstractmethod
    def perimeter(self) -> NDArray:
        """Perimeter for each dot"""

    @abstractmethod
    def round_values(self, decimals, int_type) -> None:
        """Round all values

        Args:
            decimals: number of decimal places
            int_type: numpy int type (default=numpy.int16)
        """

    @ abstractmethod
    def iter(self, indices: Optional[IntOVector] = None) -> Iterator:
        """iterate over all or a part of the objects

        Parameters
        ----------
        indices: int or interable of integer

        Notes
        -----
        To iterate all object you might all use the class iterator __iter__:
        >>> for obj in my_array:
        >>>    print(obj)
        """

    @ abstractmethod
    def find(self, diameter: Optional[float] = None,
             attribute: Any = None) -> List[int]:
        """Search for an object

        Args:
            diameter: diameter to search for (optional)
            attribute: attribute to search for (optional)

        Returns:
            indices of found objects
        """

    @ abstractmethod
    def to_dict(self) -> dict:
        """dict representation of the object area"""

    @ classmethod
    @ abstractmethod
    def from_dict(cls, the_dict: Dict[str, Any]):
        """read Dot collection from dict"""

    # generic methods
    def join(self, object_array) -> None:
        """joins with another objects list """
        self.add(object_array.iter())

    def get(self, indices: Optional[IntOVector] = None) -> List:
        """Returns list of rectangles/pictures)"""
        return list(self.iter(indices=indices))

    @abstractmethod
    def copy(self, indices, deep_copy) -> Any:
        """A (deep) copy of the nsn stimulus.

        It allows to copy a subset of dot only.

        Args:
            indices: Arrays with the indices to be be copied. If None, all
                object will be copies

        Returns:
            a copy of the array
        """
