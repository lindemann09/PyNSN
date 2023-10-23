"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from collections import OrderedDict
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union
import shapely
import numpy as np
from numpy.typing import NDArray

from .._lib.geometry import round2
from .shapes import Dot, Picture, Rectangle

IntOVector = Union[int, Sequence[int], NDArray[np.int_]]
ShapeType = Union[Dot, Rectangle, Picture]


class ShapeArray(object):
    """Numpy Optimizes Representation of Array of shapes"""

    def __init__(self) -> None:
        self._xy = np.empty((0, 2), dtype=np.float64)
        self._dot_diameter = np.empty(0, dtype=np.float64)
        self._rect_sizes = np.empty((0, 2), dtype=np.float64)
        self._polygons = np.empty(0, dtype=object)
        self._attributes = np.empty(0, dtype=object)
        self._dot_ids = np.empty(0, dtype=bool)
        self._rect_ids = np.empty(0, dtype=bool)

    @property
    def xy(self) -> NDArray:
        return self._xy

    @property
    def dot_diameter(self) -> NDArray:
        return self._dot_diameter

    @property
    def rect_sizes(self) -> NDArray:
        return self._rect_sizes

    @property
    def polygons(self) -> NDArray[shapely.Polygon]:
        return self._polygons

    @property
    def attributes(self) -> NDArray:
        return self._attributes

    @property
    def dot_ids(self) -> NDArray[np.bool_]:
        return self._dot_ids

    @property
    def rect_ids(self) -> NDArray[np.bool_]:
        return self._rect_ids

    def _find_shape_ids(self):
        self._rect_ids = np.logical_not(np.isnan(self._rect_sizes[:, 0]))
        self._dot_ids = np.logical_not(np.isnan(self._dot_diameter))

    def add(self, shapes: Union[ShapeType, Tuple, Sequence, ShapeArray]):
        if isinstance(shapes, ShapeType):
            if isinstance(shapes, Dot):
                dia = shapes.diameter
                size = np.full((1, 2), np.nan)
            elif isinstance(shapes, Rectangle):
                dia = np.nan
                size = np.atleast_2d(shapes.size)
            else:
                raise TypeError(f"Unknown ShapeType. type={type(shapes)}")

            poly = shapes.polygon
            self._polygons = np.append(self._polygons, poly)
            self._xy = np.append(self._xy, np.atleast_2d(shapes.xy), axis=0)
            self._attributes = np.append(self._attributes, shapes.attribute)
            self._rect_sizes = np.append(self._rect_sizes, size, axis=0)
            self._dot_diameter = np.append(self._dot_diameter, dia)
            self._find_shape_ids()

        elif isinstance(shapes, (list, tuple)):
            for x in shapes:
                self.add(x)

        elif isinstance(shapes, ShapeArray):
            self.add(shapes.get_list())

        else:
            raise TypeError(
                f"Can't add '{type(shapes)}', that's not a ShapeType or list of ShapeTypes."
            )

    def replace(self, index: int, shape: ShapeType):
        if isinstance(shape, Dot):
            dia = self._dot_diameter
            size = [np.nan, np.nan]
        elif isinstance(shape, Rectangle):
            dia = np.nan
            size = shape.size
        else:
            raise TypeError(f"Unknown ShapeType. type={type(shape)}")

        poly = shape.polygon
        self._polygons[index] = poly
        self._xy[index, :] = shape.xy
        self._attributes[index] = shape.attribute
        self._dot_diameter[index] = dia
        self._rect_sizes[index, :] = size
        self._find_shape_ids()


    def delete(self, index: IntOVector) -> None:
        self._polygons = np.delete(self._polygons, index)
        self._xy = np.delete(self._xy, index, axis=0)
        self._attributes = np.delete(self._attributes, index)
        self._dot_diameter = np.delete(self._dot_diameter, index)
        self._rect_sizes = np.delete(self._rect_sizes, index, axis=0)
        self._find_shape_ids()


    def clear(self):
        """ """
        self._polygons = np.empty(0, dtype=object)
        self._xy = np.empty((0, 2), dtype=np.float64)
        self._attributes = np.empty(0, dtype=object)
        self._dot_diameter = np.empty(0, dtype=np.float64)
        self._rect_sizes = np.empty((0, 2), dtype=np.float64)
        self._find_shape_ids()


    def round_values(self, decimals: int = 0, int_type: type = np.int64,
                     rebuild_polygons=True) -> None:
        """rounds all values"""

        if decimals is None:
            return
        self._xy = round2(self._xy, decimals=decimals, int_type=int_type)
        self._dot_diameter = round2(self._dot_diameter, decimals=decimals,
                                    int_type=int_type)
        self._rect_sizes = round2(self._rect_sizes, decimals=decimals,
                                  int_type=int_type)

        if rebuild_polygons:
            for i, shape in enumerate(self.get_list()):
                shape.delete_polygon()
                self._polygons[i] = shape.polygon

    def get(self, index: int) -> Union[Dot, Rectangle]:
        """Returns  selected object"""
        if np.isnan(self._dot_diameter[index]):
            rtn = Rectangle(
                xy=self._xy[index, :],
                size=self._rect_sizes[index, :],
                attribute=self._attributes[index])
        else:
            rtn = Dot(
                xy=self._xy[index, :],
                diameter=self._dot_diameter[index],
                attribute=self._attributes[index])

        rtn.set_polygon(self._polygons[index])
        return rtn

    def get_list(
        self, index: Optional[IntOVector] = None
    ) -> List[Union[Dot, Rectangle]]:
        """Returns list with selected objects
        """
        if index is None:
            ids = range(self.n_objects)
        elif isinstance(index, int):
            ids = (index,)
        else:
            ids = list(index)

        return [self.get(x) for x in ids]

    @property
    def n_objects(self) -> int:
        """number of shapes"""
        return len(self._attributes)

    def to_dict(self) -> OrderedDict:
        """dict representation of the object array

        Notes. The can be used to create Pandas dataframe or Arrow Tables

        Examples
        --------
        >>> df_dict = stimulus.dataframe_dict()
        >>> df = pandas.DataFame(df_dict) # Pandas dataframe

        >>> table = pyarrow.Table.from_pydict(df_dict) # Arrow Table
        """

        d = OrderedDict()
        d.update({"x": self._xy[:, 0].tolist(),
                  "y": self._xy[:, 1].tolist(),
                  "diameter": self._dot_diameter.tolist(),
                  "width": self._rect_sizes[:, 0].tolist(),
                  "height": self._rect_sizes[:, 1].tolist(),
                  "attributes": self._attributes.tolist()})
        return d

    def csv(self, skip_columns: Optional[Sequence[str]] = None) -> str:
        """Comma-separated table representation the nsn stimulus

        Args:
            variable_names: if True first line include the variable names
            hash_column: if True hash will be included as first column
            skip_columns: list of column names that should not be exported
        Returns:
            CSV representation

        """

        d = self.to_dict()
        rtn = ",".join(d.keys()) + "\n"
        keys = list(d.keys())
        if skip_columns is not None:
            for sc in skip_columns:
                try:
                    keys.remove(sc)
                except ValueError:
                    pass

        for i in range(self.n_objects):
            row = ""
            for k in keys:
                row += f"{d[k][i]},"
            row = row[:-1] + "\n"  # replace comma
            rtn += row
        return rtn

    # def iter(self) -> Iterator:
    #     """iterate over all objects

    #     Parameters
    #     ----------
    #     indices: int or iterable of integer

    #     Notes
    #     -----
    #     To iterate all object you might all use the class iterator __iter__:
    #     >>> for obj in my_array:
    #     >>>    print(obj)
    #     """
    #     pass


def from_dict(the_dict: Dict[str, Any]) -> ShapeArray:
    """read shape array from dict"""
    ##
    pass  # FIXME
