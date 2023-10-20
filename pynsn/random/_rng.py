"""NOTE:
import always the entire module and call `rng.generator` to ensure that you
access of the newly initialized random generator, after calling `init_random_generator`
"""

from typing import Optional, Sequence, Union

import numpy as np
import shapely
from numpy.typing import NDArray

generator = np.random.default_rng()


def init_random_generator(seed: Union[int, NDArray, Sequence, None] = None):
    """Init random generator and set random seed (optional)

    Parameters
    ----------
    seed: seed value
        must be none, int or array_like[ints]

    Notes
    -----
    see documentation of `numpy.random.default_rng()` of Python standard library
    """

    # pylint:disable = W0603
    global generator
    generator = np.random.default_rng(seed=seed)
    if seed is not None:
        print(f"PyNSN seed: {seed}")


class BrownianMotion(object):
    """Random walk of shapley.Points or Polygons"""

    def __init__(self,
                 point_or_polygon: Union[shapely.Point, shapely.Polygon],
                 walk_area: Optional[shapely.Polygon] = None,
                 delta: float = 2,
                 bounce_boarder=True):
        """performs brownian motion (random walk) of an position of shapely.polygon

        walk area can be limited.

        Parameters
        ----------

        point_or_polygon: coordinate or shapley.polygon
            Point or Polygon to move

        bounce_boarder: bool
            if true, random walk bounces back at walk_area boarder, otherwise walk
            will be continued in the initial start position

        Notes
        -----
        see Brownian motion https://en.wikipedia.org/wiki/Brownian_motion
        """

        if isinstance(walk_area, shapely.Polygon):
            shapely.prepare(walk_area)
            if not walk_area.contains_properly(point_or_polygon):
                raise ValueError("start shape is outside walk area")
            self.area = walk_area
        else:
            self.area = None  # unbounded walk

        self._geom = point_or_polygon
        self._start_geom = self._geom

        self.scale = delta ** 2
        self.bounce = bounce_boarder
        self.counter: int = 0
        self._position = np.zeros(2)

    @property
    def walk(self) -> NDArray:
        """current position of the walk in space relative to the start

        returns numpy array with xy coordinates
        """
        return self._position

    @property
    def polygon(self) -> Union[shapely.Point, shapely.Polygon]:
        """current polygon in the random walk"""
        return self._geom

    @property
    def centriod(self) -> NDArray:
        """current centroid (position)

        returns numpy array with xy coordinates
        """
        return shapely.get_coordinates(self._geom.centroid)

    def next(self, dt=1) -> None:
        """do the next step in the random walk"""
        self.counter += 1
        while True:
            step = generator.normal(loc=0, scale=self.scale * dt, size=2)
            new_geom = shapely.transform(self._geom, lambda x: x + step)
            if not isinstance(self.area, shapely.Polygon) or \
                    self.area.contains_properly(new_geom):
                # unbounded area or inside area -> good step
                self._geom = new_geom
                self._position = self._position + step
                return
            elif not self.bounce:
                # start from start position
                self._geom = self._start_geom
                self._position = np.zeros(2)
                return
