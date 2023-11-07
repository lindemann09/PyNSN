"""

"""
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
import warnings
from hashlib import md5
from typing import Optional, Tuple, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .. import constants
from .._misc import key_value_format
from .._shapes import Dot, PolygonShape, Rectangle, Ellipse, ShapeType, Point2D
from ..types import NoSolutionError
from .properties import ArrayProperties
from .shape_array import ShapeArray
from ..random import generator


class NSNStimulus(ShapeArray):
    """Non-Symbolic Number Stimulus

    NSN-Stimulus are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 target_area: Union[Dot, Rectangle, Ellipse],
                 min_distance: int = 2,
                 min_distance_target_area: int = 2
                 ) -> None:
        super().__init__()
        if not isinstance(target_area, (Dot, Rectangle, Ellipse)):
            raise RuntimeError(f"Target area can not be a {type(target_area)}")
        if (target_area.height < 1 or target_area.width < 1):
            raise RuntimeError(
                f"Target area is too small. size={target_area.size}")

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

    def fix_overlap(self,
                    inside_convex_hull: bool = False,
                    minimal_replacing: bool = True,
                    sort_before: bool = True) -> bool:  # FIXME manipulation objects
        """move an selected object that overlaps to an free position in the
        neighbourhood.

        minimal_replacing: try to find a new random position is a neighbourhood,
            otherwise overlapping object will be randomly replaced anywere in the
            search area
        returns True if position has been changed

        raise exception if not found
        occupied space: see generator generate
        """
        if sort_before:
            self.sort_by_excentricity()

        if inside_convex_hull:
            area = PolygonShape(self.convex_hull.polygon)
            area_ring = area.polygon.exterior
            shapely.prepare(area_ring)
            raise NotImplementedError("Needs testing")  # FIXME
        else:
            area = self.target_area
            area_ring = self._area_ring

        changes = False

        cnt = 0
        while cnt < 20:
            resp = np.empty(0, dtype=int)
            for x in range(self.n_objects):
                r = self._fix_overlap(index=x, min_distance=self.min_distance,
                                      target_area=area, target_area_ring=area_ring,
                                      minimal_replacing=minimal_replacing,
                                      min_distance_target_area=self.min_distance_target_area)
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
        return super().shape_overlaps(shape, min_distance)

    def matrix_overlaps(self, min_distance: Optional[float] = None) -> NDArray:
        if min_distance is None:
            min_distance = self.min_distance
        return self.matrix_dwithin(distance=min_distance)

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

        while n > 0:
            try:
                new_object = self.random_free_position(
                    shape=ref_object.copy(),
                    ignore_overlaps=ignore_overlaps,
                    inside_convex_hull=inside_convex_hull)
            except NoSolutionError as err:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={self.n_objects}") from err

            self.add(new_object)
            n = n - 1

    def random_free_position(self,
                             shape: ShapeType,
                             ignore_overlaps: bool = False,
                             inside_convex_hull: bool = False) -> ShapeType:
        """moves the object to a random free position

        raise exception if not found
        """
        if not isinstance(shape, (ShapeType)):
            raise NotImplementedError("Not implemented for "
                                      f"{type(shape).__name__}")

        if inside_convex_hull:
            area = PolygonShape(self.convex_hull.polygon)
            area_ring = area.polygon.exterior
            shapely.prepare(area_ring)
        else:
            area = self.target_area
            area_ring = self._area_ring

        return self._random_free_position(shape=shape,
                                          min_distance=self.min_distance,
                                          min_distance_target_area=self.min_distance_target_area,
                                          ignore_overlaps=ignore_overlaps,
                                          target_area=area,
                                          target_area_ring=area_ring)
