"""

"""
from __future__ import annotations
from copy import deepcopy

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
import warnings
from hashlib import md5
from typing import Tuple, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .. import constants
from .._misc import key_value_format
from .._shapes import Dot, PolygonShape, Rectangle, ShapeType, Point2D
from ..random._rng import WalkAround, generator
from ..types import NoSolutionError
from .properties import ArrayProperties
from .shape_array import ShapeArray


class NSNStimulus(ShapeArray):
    """Non-Symbolic Number Stimulus

    NSN-Stimulus are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 target_area: Union[Dot, Rectangle],
                 min_distance: int = 2,
                 min_distance_target_area: int = 2
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

        self.target_area = target_area
        if tuple(target_area.xy) != (0, 0):
            warnings.warn("TargetArea does not use shape position. "
                          "Shape Position will be set to (0, 0).",
                          UserWarning)
            self.target_area.xy = (0, 0)

        self._area_ring = shapely.get_exterior_ring(self.target_area.polygon)
        shapely.prepare(self.target_area.polygon)
        shapely.prepare(self._area_ring)

        self.min_distance = min_distance
        self.min_distance_target_area = min_distance_target_area
        self._properties = ArrayProperties(self)

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

    def properties_txt(self, with_hash: bool = False, extended_format: bool = False) -> str:
        prop = self.properties
        if extended_format:
            if with_hash:
                rtn = f"- Hash {self.hash()}\n "
            else:
                rtn = ""
            first = True
            for k, v in prop.to_dict().items():
                if first and len(rtn) == 0:
                    rtn = "- "
                    first = False
                else:
                    rtn += " "
                rtn += key_value_format(k, v) + "\n "

        else:
            if with_hash:
                rtn = "HASH: {} ".format(self.hash())
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

    # def random_position(self) -> NDArray:
    #     """random position inside the target area"""
    #     b = self.target_area.polygon.bounds  # l, b, r, t
    #     bounds_size = np.array([b[2]-b[0], b[3]-b[1]])

    #     while True:
    #         pos = rng.generator.random(size=2) * bounds_size - (bounds_size/2)
    #         if self.target_area.polygon.contains_properly(shapely.Point(pos)):
    #             return pos



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

    def fix_overlap(self, object_index: int,
                    inside_convex_hull: bool = False) -> bool:  # FIXME manipulation objects
        """move an selected object that overlaps to an free position in the
        neighbourhood.

        returns True if position has been changed

        raise exception if not found
        occupied space: see generator generate
        """

        target = self.pop(object_index)

        if inside_convex_hull:
            search_area = self.convex_hull.polygon
        else:
            search_area = self.target_area.polygon

        random_walk = WalkAround(target.polygon,
                                 walk_area=search_area,
                                 delta=2)

        while True:
            # polygons list is prepared
            overlaps = shapely.dwithin(
                self.polygons, random_walk.polygon, distance=0)
            if not np.any(overlaps):
                break  # good position
            if random_walk.counter > constants.MAX_ITERATIONS:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={self.n_objects}")
            random_walk.next()

        changes = random_walk.counter > 0
        if changes:
            target.xy = np.array(target.xy) + random_walk.walk

        self.add(target)
        return changes

    @property
    def has_overlaps(self) -> bool:
        """Returns True for two or more elements overlap (i.e. taking
        into account the minimum distance).
        """
        arr = deepcopy(self)
        while arr.n_objects>0:
            shape = arr.pop()
            y = arr.dwithin(shape=shape, distance=self.min_distance)
            if np.any(y!=0):
                return True

        return False

    def shape_overlaps(self, shape: Union[Point2D, ShapeType]) -> NDArray[np.int_]:
        """Returns True for all elements that overlap  with the particular shape
        (i.e. taking into account the minimum distance).
        """
        return self.dwithin(shape, distance=self.min_distance)

    def matrix_overlaps(self) -> NDArray:
        return self.matrix_dwithin(distance=self.min_distance)

    def inside_target_area(self, shape: Union[Point2D, ShapeType]) -> bool:
        """Returns True if shape is inside target area.
        """
        return shape.is_inside(shape=self.target_area,
                               shape_exterior_ring=self._area_ring,
                               min_dist_boarder=self.min_distance_target_area)

    def add_somewhere(self,
                      ref_object: ShapeType,
                      n: int = 1,
                      ignore_overlaps: bool = False,
                      inside_convex_hull: bool = False):
        # TODO could be improved by return a Shape with created polygon
        # if inside_convex_hull:
        #     search_area = self.properties.convex_hull.polygon
        #     search_area_ring = search_area.exterior
        #     shapely.prepare(search_area)
        #     shapely.prepare(search_area_ring)
        # else:
        #     search_area = self.target_area.polygon
        #     search_area_ring = self._area_ring
        # FIXME inside convex hull not yet implemented, Problem: ch must be centred

        if not isinstance(ref_object, (ShapeType)):
            raise NotImplementedError("Not implemented for "
                                      f"{type(ref_object).__name__}")
        while n > 0:
            try:
                new_object = self._move_to_free_position(shape=ref_object,
                                                         ignore_overlaps=ignore_overlaps,
                                                         inside_convex_hull=inside_convex_hull)
            except NoSolutionError as err:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={self.n_objects}") from err

            self.add(new_object)
            n = n - 1

    def random_free_position(self,
                             ref_object: ShapeType,
                             ignore_overlaps: bool = False,
                             inside_convex_hull: bool = False) -> Tuple[float, float]:
        """returns random free position for this object or polygon

        raise exception if not found
        occupied space: see generator generate
        """
        if not isinstance(ref_object, (ShapeType)):
            raise NotImplementedError("Not implemented for "
                                      f"{type(ref_object).__name__}")
        new_object = self._move_to_free_position(shape=ref_object,
                                                 ignore_overlaps=ignore_overlaps,
                                                 inside_convex_hull=inside_convex_hull)
        return new_object.xy

    def _move_to_free_position(self,
                               shape: ShapeType,
                               ignore_overlaps: bool = False,
                               inside_convex_hull: bool = False) -> ShapeType:
        """returns shape moved to a new free position in the target area
        """

        if inside_convex_hull:
            area = PolygonShape(self.convex_hull.polygon)
            area_ring = area.polygon.exterior
            shapely.prepare(area_ring)
            raise NotImplementedError("Needs testing")  # FIXME
        else:
            area = self.target_area
            area_ring = self._area_ring

        b = area.polygon.bounds  # l, b, r, t
        search_bounds_size = np.array([b[2]-b[0], b[3]-b[1]])

        candidate = shape.copy()
        cnt = 0
        while True:
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError(
                    "Can't find a free position for this polygon")
            cnt += 1
            # propose a random position
            candidate.xy = generator.random(size=2) * search_bounds_size \
                - (search_bounds_size/2)

            if not candidate.is_inside(shape=area,
                                       shape_exterior_ring=area_ring,
                                       min_dist_boarder=self.min_distance_target_area):
                continue

            if ignore_overlaps:
                return candidate
            else:
                # find overlaps
                overlaps = self.dwithin(candidate, distance=self.min_distance)
                if np.all(overlaps==0):
                    return candidate
