"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
import warnings
from hashlib import md5
from typing import Optional, Sequence, Tuple, Union
from numpy.typing import NDArray
import numpy as np
import shapely

from .. import constants
from .._lib.exceptions import NoSolutionError
from .._lib.misc import key_value_format
from ..random._rng import BrownianMotion, generator
from .properties import ArrayProperties
from .shape_array import IntOVector, ShapeArray, ShapeType
from .shapes import Dot, Picture, Rectangle


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

    def copy(self) -> NSNStimulus:
        """A copy of the nsn stimulus.

        Returns:
            a copy of the array
        """
        rtn = NSNStimulus(target_area=self.target_area.copy(),
                          min_distance=self.min_distance)
        rtn.add(self.get_list())  # copy all objects
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
    #         if self.target_area.polygon.contains_properly(shapely.Point(pos)):
    #             return pos

    def add(self, shapes: Union[ShapeType, Tuple, Sequence, ShapeArray]):
        """add shape a the defined position

        Note
        ----
        see `add_somewhere` for adding shape at a random position
        """
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

    def add_somewhere(self,
                      ref_object: ShapeType,
                      n: int = 1,
                      ignore_overlaps: bool = False,
                      inside_convex_hull: bool = False):

        if not isinstance(ref_object, (Dot, Rectangle, Picture)):
            raise NotImplementedError("Not implemented for "
                                      f"{type(ref_object).__name__}")
        ctr_object = ref_object.copy()
        ctr_object.xy = (0, 0)  # ensure centred position
        while n > 0:
            new_object = ctr_object.copy()
            new_object.xy = self._random_free_position(
                centred_polygon=ctr_object.polygon,
                ignore_overlaps=ignore_overlaps,
                inside_convex_hull=inside_convex_hull)
            self.add(new_object)
            n = n - 1

    def _random_free_position(self,
                              centred_polygon: shapely.Polygon,
                              ignore_overlaps: bool = False,
                              inside_convex_hull: bool = False) -> NDArray:
        """returns random free position for this object or polygon

        raise exception if not found
        occupied space: see generator generate
        """
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
        if inside_convex_hull:
            raise NotImplementedError()

        if ignore_overlaps:
            other_polygons = None
        else:
            other_polygons = self.polygons

        try:
            xy = move_to_free_position(polygon=centred_polygon,
                                       target_area=self.target_area.polygon,
                                       other_polygons=other_polygons,
                                       min_distance=self.min_distance,
                                       min_dist_area_boarder=self.min_distance_target_area,
                                       target_area_ring=self._area_ring)
        except NoSolutionError as err:
            raise NoSolutionError("Can't find a free position: "
                                  + f"Current n={self.n_objects}") from err

        return xy

    def random_free_position(self,
                             ref_object: ShapeType,
                             ignore_overlaps: bool = False,
                             inside_convex_hull: bool = False) -> NDArray:
        """returns random free position for this object or polygon

        raise exception if not found
        occupied space: see generator generate
        """
        ctr_object = ref_object.copy()
        ctr_object.xy = (0, 0)  # ensure centred position
        return self._random_free_position(centred_polygon=ctr_object.polygon,
                                          ignore_overlaps=ignore_overlaps,
                                          inside_convex_hull=inside_convex_hull)

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
            search_area = self.properties.convex_hull.polygon
        else:
            search_area = self.target_area.polygon

        random_walk = BrownianMotion(
            target.make_polygon(buffer=self.min_distance),
            walk_area=search_area, delta=2)

        while True:
            # polygons list is prepared
            overlaps = shapely.intersects(self.polygons, random_walk.polygon)
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

    def get_overlaps(self, polygon: shapely.Polygon):

        return shapely.dwithin(self.polygons, polygon, distance=self.min_distance)


def move_to_free_position(polygon: shapely.Polygon,
                          target_area: shapely.Polygon,
                          other_polygons: Optional[NDArray[shapely.Polygon]] = None,
                          min_distance: int = 0,
                          min_dist_area_boarder: Optional[int] = None,
                          target_area_ring: Optional[shapely.Polygon] = None
                          ) -> NDArray:
    """returns the required movement (x, y) of the polygon to a randomly chosen
    position in the target area.

    target_ring: the exterior of target_area
        it can be explicitly specified to avoid generation of exterior polygon and
        thus improve performance with reoccurring function calls

    Note:


     search performance is better if `target_area`, `target_area_ring` and
     `other_polygons` are prepared (see `shapely.prepare()`)
    """

    if not shapely.contains_properly(target_area, polygon):
        raise RuntimeError("Polygon has to been inside target_area")

    if min_dist_area_boarder is None:
        min_dist_area_boarder = min_distance
    if target_area_ring is None:
        target_area_ring = target_area.exterior

    b = target_area.bounds  # l, b, r, t
    search_bounds_size = np.array([b[2]-b[0], b[3]-b[1]])
    movement = np.zeros(2)
    cnt = 0
    while True:
        if cnt > constants.MAX_ITERATIONS:
            raise NoSolutionError(
                "Can't find a free position for this polygon")
        cnt += 1

        # propose a random position
        movement = generator.random(size=2) * search_bounds_size \
            - (search_bounds_size/2)
        candidate_polygon = shapely.transform(polygon,
                                              lambda x: x + movement)
        # target is inside if covered and not too close to target_area_ring
        if target_area.contains_properly(candidate_polygon) and not \
            shapely.dwithin(target_area_ring, candidate_polygon,
                            distance=min_dist_area_boarder):

            if other_polygons is None:
                return movement
            else:
                # find overlaps
                overlaps = shapely.dwithin(
                    candidate_polygon, other_polygons, distance=min_distance)
                if not np.any(overlaps):
                    return movement
