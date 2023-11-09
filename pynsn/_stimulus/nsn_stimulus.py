"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
from hashlib import md5
from typing import Optional, Union

import numpy as np
from numpy.typing import NDArray

from .._shapes import Dot, PolygonShape, Rectangle, Ellipse, ShapeType, Point2D
from ..types import NoSolutionError
from .properties import ArrayProperties
from .shape_array import ShapeArray
from .target_area import TargetArea
from .stimulus_colours import StimulusColours
from .. import defaults
from ..random import MultiVarDistributionType


class NSNStimulus(ShapeArray):
    """Non-Symbolic Number Stimulus

    NSN-Stimulus are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 target_area_shape: Union[Dot, Rectangle, Ellipse, PolygonShape],
                 min_distance: int = defaults.MIN_DISTANCE,
                 min_distance_target_area: int = defaults.MIN_DISTANCE,
                 random_distribution: Optional[MultiVarDistributionType] = None
                 ) -> None:

        super().__init__()
        self._target_area = TargetArea(shape=target_area_shape,
                                       min_distance_boarder=min_distance_target_area,
                                       random_distribution=random_distribution)
        self.min_distance = min_distance
        self._properties = ArrayProperties(self)
        self._colours = StimulusColours(target_area=self._target_area.colour)

    @property
    def target_area(self) -> TargetArea:
        """the target area of the stimulus"""
        return self._target_area

    @property
    def colours(self) -> StimulusColours:
        """the colours of the stimulus"""
        return self._colours

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

    def properties_txt(self, with_hash: bool = False, short_format: bool = False) -> str:
        if with_hash:
            if not short_format:
                rtn = f"- Hash {self.hash()}\n "
            else:
                rtn = "HASH: {} ".format(self.hash())
        else:
            rtn = ""

        return rtn + self._properties.totext(short_format)

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
        rtn.update(self._attributes.tobytes())
        return rtn.hexdigest()

    def fix_overlap(self,
                    inside_convex_hull: bool = False,
                    minimal_replacing: bool = True,
                    sort_before: bool = True,
                    max_iterations: Optional[int] = None) -> bool:  # FIXME manipulation objects
        """move an selected object that overlaps to an free position in the
        neighbourhood.

        minimal_replacing: try to find a new random position is a neighbourhood,
            otherwise overlapping object will be randomly replaced anywere in the
            search area
        returns True if position has been changed

        raise exception if not found
        occupied space: see generator generate
        """
        if max_iterations is None:
            max_iterations = defaults.MAX_ITERATIONS

        if sort_before:
            self.sort_by_excentricity()

        if inside_convex_hull:
            area = TargetArea(
                shape=PolygonShape(self.convex_hull.polygon),
                min_distance_boarder=self._target_area.min_distance_boarder)
        else:
            area = self._target_area

        changes = False
        cnt = 0
        while cnt < 20:
            resp = np.empty(0, dtype=int)
            for x in range(self.n_objects):
                r = self._fix_overlap(index=x,
                                      min_distance=self.min_distance,
                                      minimal_replacing=minimal_replacing,
                                      target_area=area,
                                      max_iterations=max_iterations)
                resp = np.append(resp, r)
            if np.any(resp == 1):
                changes = True
            if not np.any(resp == -1):  # solution found?
                return changes
            cnt += 1

        raise NoSolutionError("Can't find a solution with no overlaps")

    def contains_overlaps(self, min_distance: Optional[float] = None) -> bool:
        """Returns True for two or more elements overlap (i.e. taking
        into account the minimum distance).
        """
        if min_distance is None:
            min_distance = self.min_distance
        return super().contains_overlaps(min_distance)

    def get_overlaps(self, index: int,
                     min_distance: Optional[float] = None) -> NDArray[np.bool_]:
        """get overlaps with other objects. Ignores overlap with oneself."""
        if min_distance is None:
            min_distance = self.min_distance
        return super().get_overlaps(index, min_distance)

    def shape_overlaps(self, shape: Union[Point2D, ShapeType],
                       min_distance: Optional[float] = None) -> NDArray[np.bool_]:
        """Returns True for all elements that overlap  with the particular shape
        (i.e. taking into account the minimum distance).
        """
        if min_distance is None:
            min_distance = self.min_distance
        return self.dwithin(shape, distance=min_distance)

    def matrix_overlaps(self, min_distance: Optional[float] = None) -> NDArray:
        if min_distance is None:
            min_distance = self.min_distance
        return self.matrix_dwithin(distance=min_distance)

    def inside_target_area(self, shape: Union[Point2D, ShapeType]) -> bool:
        """Returns True if shape is inside target area.
        """
        return self._target_area.is_object_inside(shape)

    def add_somewhere(self,
                      ref_object: ShapeType,
                      n: int = 1,
                      ignore_overlaps: bool = False,
                      inside_convex_hull: bool = False,
                      max_iterations: Optional[int] = None):

        while n > 0:
            try:
                new_object = self.random_free_position(
                    shape=ref_object.copy(),
                    ignore_overlaps=ignore_overlaps,
                    inside_convex_hull=inside_convex_hull,
                    max_iterations=max_iterations)
            except NoSolutionError as err:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={self.n_objects}") from err

            self.add(new_object)
            n = n - 1

    def random_free_position(self,
                             shape: ShapeType,
                             ignore_overlaps: bool = False,
                             inside_convex_hull: bool = False,
                             max_iterations: Optional[int] = None) -> ShapeType:
        """moves the object to a random free position

        raise exception if not found
        """
        if not isinstance(shape, (ShapeType)):
            raise NotImplementedError("Not implemented for "
                                      f"{type(shape).__name__}")
        if max_iterations is None:
            max_iterations = defaults.MAX_ITERATIONS

        if inside_convex_hull:
            area = TargetArea(
                shape=PolygonShape(self.convex_hull.polygon),
                min_distance_boarder=self._target_area.min_distance_boarder)
        else:
            area = self._target_area

        return self._random_free_position(shape=shape,
                                          min_distance=self.min_distance,
                                          ignore_overlaps=ignore_overlaps,
                                          target_area=area,
                                          max_iterations=max_iterations)
