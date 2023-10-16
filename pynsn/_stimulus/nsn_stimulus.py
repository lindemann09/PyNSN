"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
from numpy.typing import NDArray

import warnings
from typing import Union
import shapely
import numpy as np

from .shapes import Dot, Rectangle, Picture
from .properties import ArrayProperties
from .shape_array import ShapeArray
from .. import _stimulus, constants
from .._lib.misc import key_value_format
from .._lib import rng
from .._lib.exceptions import NoSolutionError


class NSNStimulus(object):
    """Non-Symbolic Number Stimulus

    NSN-Stimulus are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(
        self,
        target_area: Union[Dot, Rectangle],
        min_dist_between_shapes: float = 2,
        min_dist_area_edge: float = 2,
    ) -> None:
        super().__init__()
        assert isinstance(target_area, (Dot, Rectangle))
        if isinstance(target_area, Dot) and target_area.diameter < 1:
            raise RuntimeError(
                f"Target area is too small. Diameter={target_area.diameter}")
        if isinstance(target_area, Rectangle) and \
                (target_area.height < 1 or target_area.width < 1):
            raise RuntimeError(
                f"Target area is too small. Size={target_area.size}")

        if tuple(target_area.xy) != (0, 0):
            warnings.warn("TargetArea does not use shape position. "
                          "Shape Position will be set to (0, 0).",
                          UserWarning)

        self.target_area = target_area
        self.target_area.xy = (0, 0)
        shapely.prepare(self.target_area.polygon)

        self.min_dist_between_objects = min_dist_between_shapes
        self.min_dist_area_edge = min_dist_area_edge
        self._shapes = ShapeArray()

    @property
    def shapes(self) -> ShapeArray:
        """ShapeArray of the stimulus

        see documentation ``ShapeArray``
        """
        return self._shapes

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
        return self._shapes.properties

    def copy(self) -> NSNStimulus:
        """A copy of the nsn stimulus.

        Returns:
            a copy of the array
        """
        rtn = NSNStimulus(
            target_area=self.target_area.copy(),
            min_dist_between_shapes=self.min_dist_area_edge,
            min_dist_area_edge=self.min_dist_area_edge
        )
        rtn.shapes.add(self._shapes)
        return rtn

    def properties_txt(self, with_hash: bool = False, extended_format: bool = False) -> str:
        prop = self.properties
        if extended_format:
            rtn = None
            for k, v in prop.to_dict().items():
                if rtn is None:
                    if with_hash:
                        rtn = f"-{k}  {v}\n "
                    else:
                        rtn = "-"
                else:
                    rtn += key_value_format(k, v) + "\n "

            if rtn is None:
                rtn = ""
        else:
            if with_hash:
                rtn = "ID: {} ".format(self._shapes.hash())
            else:
                rtn = ""
            rtn += (f"N: {prop.numerosity}, "
                    + f"TSA: {int(prop.total_surface_area)}, "
                    + f"ISA: {int(prop.average_surface_area)}, "
                    + f"FA: {int(prop.field_area)}, "
                    + f"SPAR: {prop.sparsity:.2f}, "
                    + f"logSIZE: {prop.log_size:.2f}, "
                    + f"logSPACE: {prop.log_spacing:.2f}, "
                    + f"COV: {prop.coverage:.2f}")

        return rtn.rstrip()

    def random_position(self) -> NDArray:
        """random position inside the target area"""
        b = self.target_area.polygon.bounds  # l, b, r, t
        bounds_size = np.array([b[2]-b[0], b[3]-b[1]])

        while True:
            pos = rng.generator.random(size=2) * bounds_size - (bounds_size/2)
            if self.target_area.polygon.covers(shapely.Point(pos)):
                return pos

    def find_position(self,
                      ref_object: Union[Dot, Rectangle, Picture],
                      in_neighbourhood: bool = False,
                      allow_overlapping: bool = False,
                      inside_convex_hull: bool = False,
                      occupied_space=None) -> _stimulus.ShapeType:
        """returns the copy of the object of at a randomly choose free position

        raise exception if not found
        occupied space: see generator generate
        """

        if not isinstance(ref_object, (Dot, Rectangle, Picture)):
            raise NotImplementedError("Not implemented for "
                                      f"{type(ref_object).__name__}")

        if occupied_space is not None and \
                not isinstance(occupied_space, NSNStimulus):
            raise TypeError(
                "Occupied_space has to be a Dot or nsn stimulus or None.")

        if inside_convex_hull:
            search_area = self.properties.convex_hull.polygon
            shapely.prepare(search_area)
        else:
            search_area = self.target_area.polygon

        b = search_area.bounds  # l, b, r, t
        search_bounds_size = np.array([b[2]-b[0], b[3]-b[1]])

        target = ref_object.copy()
        # random position assume in contrast to random walk neighbourhood that pos is (0,0)
        if not in_neighbourhood:
            target.xy = (0, 0)
        # takes into account minimum distance for all comparisons
        target_polygon = target.make_polygon(buffer=constants.DEFAULT_MIN_DIST)
        if in_neighbourhood:
            # start random walk
            random_walk = rng.BrownianMotion(target_polygon,
                                             walk_area=search_area, delta=2)
        else:
            random_walk = None
        candidate_pos = np.zeros(2)
        cnt = 0
        while True:
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={self.shapes.n_objects}")
            cnt += 1

            if random_walk is None:
                # propose a random position
                candidate_pos = rng.generator.random(size=2) * search_bounds_size \
                    - (search_bounds_size/2)
                candidate_polygon = shapely.transform(target_polygon,
                                                      lambda x: x + candidate_pos)
                if not search_area.covers(candidate_polygon):
                    #  target not inside target area -> new proposal
                    continue
            else:
                # next random walk
                random_walk.next()
                candidate_polygon = random_walk.polygon

            if not allow_overlapping:
                # find overlaps
                for p in self.shapes.polygons:
                    if p.overlaps(candidate_polygon):  # polygons p are prepared
                        continue  # try another position
                if occupied_space is not None:
                    # check for overlap occupied space
                    raise NotImplementedError("Not yet implemented")  # FIXME

            break

        if random_walk is not None:
            target.xy = np.array(target_polygon.xy) + random_walk.walk
        else:
            target.xy = candidate_pos
        return target
