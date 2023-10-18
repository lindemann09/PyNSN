"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
import warnings
from hashlib import md5
from typing import Sequence, Tuple, Union

import numpy as np
import shapely

from .. import constants
from .._lib.exceptions import NoSolutionError
from .._lib.misc import key_value_format
from ..random._rng import BrownianMotion
from .properties import ArrayProperties
from .shape_array import IntOVector, ShapeArray, ShapeType
from .shapes import Dot, Picture, Rectangle


class NSNStimulus(ShapeArray):
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
        rtn.add(self.get_list())
        return rtn

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
    #         if self.target_area.polygon.covers(shapely.Point(pos)):
    #             return pos

    def add(self, shapes: Union[ShapeType, Tuple, Sequence, ShapeArray]):
        super().add(shapes)
        self._properties.reset_convex_hull()

    def replace(self, index: int, shape: ShapeType):
        super().replace(index, shape)
        self._properties.reset_convex_hull()

    def delete(self, index: IntOVector) -> None:
        super().delete(index)
        self._properties.reset_convex_hull()

    def pop(self, index: int) -> Union[Dot, Rectangle]:
        """Remove and return item at index"""
        rtn = self.get(index)
        self.delete(index)
        self._properties.reset_convex_hull()
        return rtn

    def clear(self):
        """ """
        super().clear()
        self._properties.reset_convex_hull()

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



    def round_values(self, decimals: int = 0, int_type: type = np.int64,
                     rebuild_polygons=True) -> None:
        """rounds all values"""
        super().round_values(decimals=decimals, int_type=int_type,
                             rebuild_polygons=rebuild_polygons)
        if rebuild_polygons:
            self._properties.reset_convex_hull()

    def find_position(self,
                      ref_object: Union[Dot, Rectangle, Picture],
                      in_neighbourhood: bool = False,
                      ignore_overlaps: bool = False,
                      inside_convex_hull: bool = False,
                      occupied_space=None) -> ShapeType:
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
        else:
            search_area = self.target_area.polygon

        shapely.prepare(search_area)
        b = search_area.bounds  # l, b, r, t
        search_bounds_size = np.array([b[2]-b[0], b[3]-b[1]])

        target = ref_object.copy()
        if not in_neighbourhood:
            # random position assume in contrast to random walk neighbourhood that pos is (0,0)
            target.xy = (0, 0)
        reference_polygon = target.make_polygon(
            buffer=constants.DEFAULT_MIN_DIST)
        if in_neighbourhood:
            # start random walk
            random_walk = BrownianMotion(reference_polygon,
                                         walk_area=search_area, delta=2)
        else:
            random_walk = None
        candidate_pos = np.zeros(2)
        cnt = 0
        while True:
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={self.n_objects}")
            cnt += 1

            if random_walk is None:
                # propose a random position
                candidate_pos = generator.random(size=2) * search_bounds_size \
                    - (search_bounds_size/2)
                candidate_polygon = shapely.transform(reference_polygon,
                                                      lambda x: x + candidate_pos)
                if not search_area.covers(candidate_polygon):
                    #  target not inside target area -> new proposal
                    continue
            else:
                # random walk
                random_walk.next()
                candidate_polygon = random_walk.polygon

            if not ignore_overlaps:
                # find overlaps
                for p in self.polygons:
                    if p.overlaps(candidate_polygon):  # polygons p are prepared
                        continue  # try another position
                if occupied_space is not None:
                    # check for overlap occupied space
                    raise NotImplementedError("Not yet implemented")  # FIXME

            break

        if random_walk is not None:
            target.xy = np.array(target.xy) + random_walk.walk
        else:
            target.xy = candidate_pos
        return target

    def fix_overlap(self, object_index: int,
                    inside_convex_hull: bool = False,
                    occupied_space=None) -> bool:
        """move an selected object that overlaps to an free position in the
        neighbourhood.

        returns True if position has been changed

        raise exception if not found
        occupied space: see generator generate
        """

        target = self.pop(object_index)

        if inside_convex_hull:
            search_area = self.properties.convex_hull.polygon
        else:
            search_area = self.target_area.polygon

        random_walk = BrownianMotion(
            target.make_polygon(buffer=constants.DEFAULT_MIN_DIST),
            walk_area=search_area, delta=2)

        while True:
            overlaps = shapely.intersects(self.polygons, random_walk.polygon) # polygons list is prepared
            if occupied_space is not None:
                # check for overlap occupied space
                raise NotImplementedError("Not yet implemented")  # FIXME
            if not np.any(overlaps):
                break # good position
            if random_walk.counter > constants.MAX_ITERATIONS:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={self.n_objects}")
            random_walk.next()

        changes = random_walk.counter>0
        if changes:
            target.xy = np.array(target.xy) + random_walk.walk

        self.add(target)
        return changes


    def get_overlaps(self, polygon: shapely.Polygon):

        rtn = shapely.distance(self.polygons, polygon)
        return rtn
        # return np.nonzero(rtn)
