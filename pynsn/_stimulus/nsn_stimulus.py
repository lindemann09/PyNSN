""" """

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import json
import warnings
from copy import deepcopy
from hashlib import md5
from pathlib import Path
from typing import Any, Dict

import numpy as np
from numpy.typing import NDArray

from .. import defaults
from .._misc import formated_json
from .._shapes import Dot, Ellipse, Point2D, PolygonShape, Rectangle, dict_to_shape
from .._shapes.abc_shapes import AbstractShape
from ..exceptions import (
    NoSolutionError,
    ShapeOutsideError,
    ShapeOutsideWarning,
    ShapeOverlapsError,
)
from ..rnd._random_shape import AbstractRndShape
from ..rnd._rng import WalkAround
from . import shape_array
from .properties import ArrayProperties
from .stimulus_colours import StimulusColours
from .target_area import TargetArea

# TODO add optional names?


class NSNStimulus(shape_array.ShapeArray):
    """Non-Symbolic Number Stimulus

    NSN-Stimulus are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(
        self,
        target_area_shape: Dot | Rectangle | Ellipse | PolygonShape,
        min_distance: int = defaults.MIN_DISTANCE,
        min_distance_target_area: int = defaults.MIN_DISTANCE,
        name: str | None = None,
    ) -> None:
        super().__init__()
        self._target_area = TargetArea(
            shape=target_area_shape, min_dist_boarder=min_distance_target_area
        )
        self.min_distance = min_distance
        self._properties = ArrayProperties(self)
        self._colours = StimulusColours(target_area=self._target_area.colour)
        self.name = name

    @property
    def target_area(self) -> TargetArea:
        """the target area of the stimulus"""
        return self._target_area

    @property
    def colours(self) -> StimulusColours:
        """the colours of the stimulus"""
        return self._colours

    @colours.setter
    def colours(self, val: StimulusColours):
        """the colours of the stimulus"""
        assert isinstance(val, StimulusColours)
        self._colours = val

    @property
    def properties(self) -> ArrayProperties:
        """Properties of the nsn stimulus.

        ``ArrayProperties`` represents and handles (fitting, scaling) visual
        properties

        * numerosity
        * average_dot_diameter/average_rectangle_size
        * total_surface_area
        * average_surface_area
        * total_perimeter
        * average_perimeter
        * field_area
        * field_area_positions
        * sparsity
        * log_spacing
        * log_size
        * coverage
        """
        return self._properties

    def properties_txt(self, with_hash: bool = True, short_format: bool = False) -> str:
        if with_hash:
            if not short_format:
                rtn = f"- Hash  {self.hash()}\n "
            else:
                rtn = "HASH: {} ".format(self.hash())
        else:
            rtn = ""

        return rtn + self._properties.to_text(short_format)[1:]

    def hash(self) -> str:
        """Hash (MD5 hash) of the array

        The hash can be used as an unique identifier of the nsn stimulus.

        Notes:
            Hashing is based on the byte representations of the positions, perimeter
            and attributes.
        """

        rtn = md5()
        # to_byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        rtn.update(self._xy.tobytes())
        try:
            rtn.update(self.properties.perimeter.tobytes())
        except AttributeError:
            pass
        return rtn.hexdigest()

    def to_dict(self, tabular: bool = True) -> dict:
        """Dict representation of the shape array"""
        if self.name is not None:
            rtn: Dict[str, Any] = {"name": self.name}

        else:
            rtn: Dict[str, Any] = {}

        rtn.update(
            {
                "hash": self.hash(),
                "target_area": self.target_area.to_dict(),
                "min_distance": self.min_distance,
                "colours": self._colours.to_dict(),
            }
        )

        if tabular:
            rtn.update({"shape_table": self.table_dict()})
        else:
            d = {"shape_array": [x.to_dict() for x in self.shapes]}
            rtn.update(d)

        return rtn

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> NSNStimulus:
        """read shape array from dict"""
        try:
            name = d["name"]
        except KeyError:
            name = None
        ta = TargetArea.from_dict(d["target_area"])
        rtn = NSNStimulus(
            target_area_shape=ta.shape,
            min_distance_target_area=ta.min_dist_boarder,
            min_distance=d["min_distance"],
            name=name,
        )
        rtn.colours = StimulusColours.from_dict(d["colours"])

        if "shape_array" in d:
            # add shapes
            for sd in d["shape_array"]:
                s = dict_to_shape(sd)
                if isinstance(s, AbstractShape):
                    rtn.shape_add(s)
        elif "shape_table" in d:
            lst = zip(
                d["shape_table"]["type"],
                d["shape_table"]["x"],
                d["shape_table"]["y"],
                d["shape_table"]["width"],
                d["shape_table"]["height"],
                d["shape_table"]["attr"],
            )
            for t, x, y, w, h, attr in lst:
                if t == Dot.shape_type():
                    s = Dot(diameter=w, xy=(x, y), attribute=attr)
                elif t == Ellipse.shape_type():
                    s = Ellipse(size=(w, h), xy=(x, y), attribute=attr)
                elif t == Rectangle.shape_type():
                    s = Rectangle(size=(w, h), xy=(x, y), attribute=attr)
                else:
                    raise NotImplementedError(f"Type: {t}")
                rtn.shape_add(s)

        return rtn

    def to_json(
        self,
        path: str | Path | None = None,
        indent: int = 2,
        tabular: bool = True,
        decimals: int | None = None,
    ) -> str:
        """JSON representation of the object

        Parameters
        ----------
        path : str | Path | None, optional
            _description_, by default None
        indent : int, optional
            _description_, by default 2
        tabular : bool, optional
            _description_, by default True
        decimals : int | None, optional
            _description_, by default None

        Returns
        -------
        str
        """
        d = self.to_dict(tabular=tabular)
        json_str = formated_json(d, indent=indent, decimals=decimals)
        if isinstance(path, (Path, str)):
            with open(path, "w", encoding=defaults.FILE_ENCODING) as fl:
                fl.write(json_str)
        return json_str

    @staticmethod
    def from_json(path: str | Path) -> NSNStimulus:
        path = Path(path)
        if not path.is_file():
            raise RuntimeError(f"Can't find {path}.")
        with open(path, "r", encoding=defaults.FILE_ENCODING) as fl:
            d = json.load(fl)
        return NSNStimulus.from_dict(d)

    def _fix_overlaps(
        self,
        index: int,
        min_distance: float,
        minimal_replacing: bool,
        target_area: TargetArea,
        max_iterations: int | None = None,
    ) -> int:
        """Move an selected object that overlaps to an free position in the
        neighbourhood.

        minimal_replacing: try to find a new random position is a neighbourhood,
            otherwise overlapping object will be randomly replaced anywhere in the
            search area


        Returns
        -------
         0: if no overlaps exist
        -1: if object overlaps, but no new position could be found
         1: if object was replaced

        occupied space: see generator generate
        """

        if not np.any(self.overlaps(index, min_distance)):
            return 0  # no overlap

        if max_iterations is None:
            max_iterations = defaults.MAX_ITERATIONS

        target = self.shapes[index]

        if minimal_replacing:
            walk = WalkAround(target.xy)
            outside_cnt = 0
            while True:
                if walk.counter > max_iterations or outside_cnt > 20:
                    return -1  # can't find a free position

                target.xy = walk.next()
                if not target_area.is_object_inside(target):
                    outside_cnt += 1
                else:
                    outside_cnt = 0
                    overlaps = self.dwithin(target, distance=min_distance)
                    overlaps[index] = False  # ignore overlap with oneself
                    if not np.any(overlaps):
                        break  # place found
        else:
            # random position anywhere
            try:
                target = self.rnd_free_pos(
                    target,
                    inside_convex_hull=False,
                    ignore_overlaps=False,
                    max_iterations=max_iterations,
                )
            except NoSolutionError:
                return -1

        self.shape_replace(index, target)
        return 1

    def fix_overlaps(
        self,
        inside_convex_hull: bool = False,
        minimal_replacing: bool = True,
        sort_before: bool = True,
        max_iterations: int | None = None,
    ) -> bool:
        """move an selected object that overlaps to an free position in the
        neighbourhood.

        minimal_replacing: try to find a new random position is a neighbourhood,
            otherwise overlapping object will be randomly replaced anywhere in the
            search area
        returns True if position has been changed

        raise exception if not found
        occupied space: see generator generate
        """

        if sort_before:
            self.sort_by_eccentricity()

        if inside_convex_hull:
            area = TargetArea(
                shape=PolygonShape(self.convex_hull.polygon),
                min_dist_boarder=self._target_area.min_dist_boarder,
            )
        else:
            area = self._target_area

        changes = False
        cnt = 0
        while cnt < 20:
            resp = np.empty(0, dtype=int)
            for x in range(len(self._shapes)):
                r = self._fix_overlaps(
                    index=x,
                    min_distance=self.min_distance,
                    minimal_replacing=minimal_replacing,
                    target_area=area,
                    max_iterations=max_iterations,
                )
                resp = np.append(resp, r)
            if np.any(resp == 1):
                changes = True
            if not np.any(resp == -1):  # solution found?
                return changes
            cnt += 1

        raise NoSolutionError("Can't find a solution with no overlaps")

    def has_overlaps(self, min_distance: float | None = None) -> bool:
        """Returns True for two or more elements overlap (i.e. taking
        into account the minimum distance).
        """
        if min_distance is None:
            min_distance = self.min_distance
        return super().has_overlaps(min_distance)

    def overlaps(
        self, index: int, min_distance: float | None = None
    ) -> NDArray[np.bool_]:
        """get overlaps with other shapes. Ignores overlap with oneself."""
        if min_distance is None:
            min_distance = self.min_distance
        return super().overlaps(index, min_distance)

    def shape_overlaps(
        self, shape: Point2D | AbstractShape, min_distance: float | None = None
    ) -> NDArray[np.bool_]:
        """Returns True for all elements that overlap with the particular shape
        (i.e. taking into account the minimum distance).
        """
        # FIXME unreliable ellipses overlapping

        if min_distance is None:
            min_distance = self.min_distance
        return self.dwithin(shape, distance=min_distance)

    def inside_target_area(self, shape: Point2D | AbstractShape) -> bool:
        """Returns True if shape is inside target area."""
        return self._target_area.is_object_inside(shape)

    def shape_add(self, shape: AbstractShape, ignore_overlaps: bool = False):
        """ "adds shape to array"""

        if not ignore_overlaps and np.any(self.shape_overlaps(shape)):
            txt = f"Shape overlaps with array. {shape}"
            if defaults.WARNINGS:
                warnings.warn(txt, ShapeOutsideWarning)
            else:
                raise ShapeOverlapsError(txt)
        if not self.target_area.is_object_inside(shape):
            txt = f"Shape outside target array. {shape}"
            if defaults.WARNINGS:
                warnings.warn(txt, ShapeOutsideWarning)
            else:
                raise ShapeOutsideError(txt)

        super().shape_add(shape)

    def __add_random_pos(
        self,
        shape: AbstractShape,
        ignore_overlaps: bool = False,
        inside_convex_hull: bool = False,
        max_iterations: int | None = None,
    ):
        """ "adds shape to random positions in the array"""

        try:
            shape = self.rnd_free_pos(
                shape=shape,
                ignore_overlaps=ignore_overlaps,
                inside_convex_hull=inside_convex_hull,
                max_iterations=max_iterations,
            )
        except NoSolutionError as err:
            raise NoSolutionError(
                "Can't find a free position: " + f"Current n={self.n_shapes}"
            ) from err

        super().shape_add(shape)

    def shape_add_random_pos(
        self,
        ref_object: AbstractShape | AbstractRndShape,
        n: int = 1,
        ignore_overlaps: bool = False,
        inside_convex_hull: bool = False,
        max_iterations: int | None = None,
    ):
        """Creates n copies of the shape(s) or n instances of the random shape(s)
        and adds them at random positions to the array (default n=1)"""

        if isinstance(ref_object, AbstractRndShape):
            for obj in ref_object.sample(n):
                self.__add_random_pos(
                    obj,
                    ignore_overlaps=ignore_overlaps,
                    inside_convex_hull=inside_convex_hull,
                    max_iterations=max_iterations,
                )
        else:
            while n > 0:
                self.__add_random_pos(
                    deepcopy(ref_object),
                    ignore_overlaps=ignore_overlaps,
                    inside_convex_hull=inside_convex_hull,
                    max_iterations=max_iterations,
                )
                n = n - 1

    def table_dict(self) -> dict:
        """Tabular representation of the array of the shapes.

        This representation can not deal with PolygonShapes. It's useful to
        create Pandas dataframe or Arrow Tables.

        Examples
        --------
        >>> df_dict = stimulus.table_dict()
        >>> df = pandas.DataFrame(df_dict) # Pandas dataframe

        >>> table = pyarrow.Table.from_pydict(df_dict) # Arrow Table
        """

        if np.any(self.shape_types() == PolygonShape.shape_type()):
            raise RuntimeError(
                "tabular shape representation can not deal with " "PolygonShapes"
            )
        d = {
            "type": self.shape_types(),
            "x": self.xy[:, 0],
            "y": self.xy[:, 1],
            "width": self.sizes[:, 0],
            "height": self.sizes[:, 1],
            "attr": [str(x) for x in self.attributes()],
        }
        return d

    def rnd_free_pos(
        self,
        shape: AbstractShape,
        ignore_overlaps: bool = False,
        inside_convex_hull: bool = False,
        max_iterations: int | None = None,
    ) -> AbstractShape:
        """moves the object to a random free position

        raise exception if not found
        """

        if not isinstance(shape, AbstractShape):
            raise NotImplementedError("Not implemented for " f"{type(shape).__name__}")
        if max_iterations is None:
            max_iterations = defaults.MAX_ITERATIONS

        if inside_convex_hull:
            area = TargetArea(
                shape=PolygonShape(self.convex_hull.polygon),
                min_dist_boarder=self.target_area.min_dist_boarder,
            )
        else:
            area = self.target_area

        cnt = 0
        while True:
            if cnt > max_iterations:
                raise NoSolutionError("Can't find a free position for this object")
            cnt += 1
            # propose a random position
            shape.xy = area.random_xy_inside_bounds()

            if not area.is_object_inside(shape):
                continue

            if ignore_overlaps:
                return shape
            else:
                # find overlaps
                overlaps = self.dwithin(shape, distance=self.min_distance)
                if not np.any(overlaps):
                    return shape
