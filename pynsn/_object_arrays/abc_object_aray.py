"""

"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from hashlib import md5
from typing import Any, Dict, Iterator, List, Optional

import numpy as np
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
    def np_append(self) -> None:
        """ """
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

    @property
    @abstractmethod
    def center_of_mass(self) -> NDArray:
        """center of mass of all objects"""

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
        indices: int or iterable of integer

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

    @abstractmethod
    def hash(self) -> str:
        """Hash (MD5 hash) of the array

        The hash can be used as an unique identifier of the nsn stimulus.

        Notes:
            Hashing is based on the byte representations of the positions, perimeter
            and attributes.
        """

    @abstractmethod
    def dataframe_dict(self, hash_column: bool = False,
                       attribute_column: bool = True) -> dict:
        """Returns dict that can be used to create Pandas dataframe or Arrow Table

        Examples
        --------
        >>> df_dict = stimulus.objects.dataframe_dict()
        >>> df = pandas.DataFame(df_dict) # Pandas dataframe

        >>> table = pyarrow.Table.from_pydict(df_dict) # Arrow Table
        """

    @abstractmethod
    def csv(self,
            variable_names: bool = True,
            hash_column: bool = False,
            attribute_column: bool = True):
        """Comma-separated table representation the nsn stimulus

        Args:
            variable_names: if True first line include the variable names
            hash_column: if True hash will be included as first column
            attribute_column: if True attributes will be included as
                    last column
        Returns:
            CSV representation

        """

    # generic functions

    def join(self, object_array) -> None:
        """joins with another objects list """
        self.add(object_array.iter())

    def get(self, indices: Optional[IntOVector] = None) -> List:
        """Returns list of rectangles/pictures)"""
        return list(self.iter(indices=indices))


def hash_array(xy: NDArray, perimeter: NDArray, attributes: NDArray) -> str:
    rtn = md5()
    # to_byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
    rtn.update(xy.tobytes())
    try:
        rtn.update(perimeter.tobytes())
    except AttributeError:
        pass
    rtn.update(attributes.tobytes())
    return rtn.hexdigest()


def make_csv(xy, size_data_dict,
             attributes=None,
             array_hash=None,
             make_variable_names=True):
    """Helper function:  makes csv for Arrays with object of different size information
    size_data_dict: keys = variable names (e.g. width, height),
                    values vector of size data
    """
    rtn = ""
    if make_variable_names:
        if array_hash:
            rtn += "hash,"
        rtn += "x,y," + ",".join(size_data_dict.keys()) + ","
        if attributes is not None:
            rtn += "attribute,"
        rtn = rtn[:-1] + "\n"  # replace comma

    size_data = np.array(list(size_data_dict.values())).T
    if attributes is None:
        attribute_vector = [None] * len(xy)  # to have something to loop
    else:
        attribute_vector = attributes

    for pos, size, attr in zip(xy, size_data, attribute_vector):
        if array_hash:
            rtn += "{0},".format(array_hash)
        rtn += "{},{},".format(pos[0], pos[1])
        for s in size:
            rtn += "{},".format(s)
        if attributes is not None:
            rtn += "{},".format(attr)
        rtn = rtn[:-1] + "\n"  # replace comma
    return rtn
