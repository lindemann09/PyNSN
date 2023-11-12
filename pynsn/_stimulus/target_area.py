"""

"""
__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

import warnings
from typing import Optional, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .._shapes import Dot, Ellipse, Point2D, PolygonShape, Rectangle, \
    ShapeType, Colour
from .. import defaults
from ..random import MultiVarDistributionType, Uniform2D
from ..types import NoSolutionError


class TargetArea(object):

    def __init__(self,
                 shape: Union[Dot, Rectangle, Ellipse, PolygonShape],
                 min_dist_boarder: int = 2
                 ) -> None:

        super().__init__()

        if not isinstance(shape, (Dot, Rectangle, Ellipse, PolygonShape)):
            raise RuntimeError(f"Target area can not be a {type(shape)}")
        if (shape.height < 1 or shape.width < 1):
            raise RuntimeError(
                f"Target area is too small. size={shape.size}")

        self._shape = shape
        if tuple(shape.xy) != (0, 0):
            warnings.warn("TargetArea does not use shape position. "
                          "Shape Position will be set to (0, 0).",
                          UserWarning)
            self._shape.xy = (0, 0)

        if self._shape.colour.value is None:
            self._shape.attribute = Colour(defaults.COLOUR_TARGET_AREA)

        self._ring = shapely.get_exterior_ring(self._shape.polygon)

        # buffer of random position samples to avoid frequent call of random distribution
        self._random_buffer = np.empty(0)

        # default random distribution is uniform
        f = np.array((-0.5, 0.5))
        self._default_distr = Uniform2D(x_minmax=shape.size[0] * f,
                                        y_minmax=shape.size[1] * f)

        self.min_dist_boarder = min_dist_boarder
        shapely.prepare(self._shape.polygon)
        shapely.prepare(self._ring)

    def is_object_inside(self, shape: Union[Point2D, ShapeType]) -> bool:
        """returns True if shape is fully inside Target area

        takes into account the minimum distance to boarder"""

        return shape.is_inside(shape=self._shape,
                               shape_exterior_ring=self._ring,
                               min_dist_boarder=self.min_dist_boarder)

    @property
    def colour(self) -> Colour:
        return self._shape.colour

    @property
    def shape(self) -> Union[Dot, Rectangle, Ellipse, PolygonShape]:
        """Shape describing the target area"""
        return self._shape

    @property
    def bound_sizes(self) -> NDArray[np.float_]:
        """width and height of bounds of the target area"""
        return self._shape.size

    def random_xy_inside_bounds(self, distr: Optional[MultiVarDistributionType] = None) -> NDArray:
        """returns a random position inside the bounds of the target area
        """
        if isinstance(distr, MultiVarDistributionType):
            return distr.sample(1)[0]
        else:
            return self._default_distr.sample(1)[0]

    def random_position(self, shape: ShapeType,
                        max_iterations: Optional[int] = None) -> ShapeType:
        """returns the object at random free position inside the target area

        takes into account the minimum distance to boarder
        """
        if max_iterations is None:
            max_iterations = defaults.MAX_ITERATIONS
        cnt = 0
        while True:
            if cnt > max_iterations:
                raise NoSolutionError(
                    "Can't find a free position for this polygon")
            cnt += 1
            # propose a random position
            shape.xy = self.random_xy_inside_bounds()

            if self.is_object_inside(shape):
                return shape

    def todict(self) -> dict:
        return {"shape": self.shape.todict(),
                "min_dist_boarder": self.min_dist_boarder}
