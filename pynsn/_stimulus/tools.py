"""
helper functions for arrays
"""
import numpy as np

from .._lib import geometry, rng


def radial_replacement_from_reference_dots(xy, ref_pos_id,
                                           neighbour_ids, replacement_size):
    """remove neighbouring position radially from reference position
    helper function, typically used for realign
    """

    # check if there is an identical position and jitter to avoid fully overlapping positions
    if np.sum(np.all(xy[neighbour_ids, ] == xy[ref_pos_id, :],
                     axis=1)) > 0:
        xy = geometry.jitter_identical_coordinates(xy)

    # relative polar positions to reference_dot
    tmp_polar = geometry.cartesian2polar(
        xy[neighbour_ids, :] - xy[ref_pos_id, :])
    tmp_polar[:, 0] = 0.000000001 + replacement_size  # determine movement size
    tmp_xy = geometry.polar2cartesian(tmp_polar)
    xy[neighbour_ids, :] = np.array([xy[neighbour_ids, 0] + tmp_xy[:, 0],
                                     xy[neighbour_ids, 1] + tmp_xy[:, 1]]).T
    return xy


def remove_overlap_from_inner_to_outer(xy, min_dist_between_objects, distance_matrix_function):
    """returns xy and boolean, if replacements were required"""
    assert callable(distance_matrix_function)

    replacement_required = False
    # from inner to outer remove overlaps
    for i in np.argsort(geometry.cartesian2polar(xy, radii_only=True)):
        dist_mtx = distance_matrix_function(between_positions=False)
        dist = dist_mtx[i, :]
        idx_overlaps = np.flatnonzero(
            dist < min_dist_between_objects).tolist()  # overlapping dot ids
        if len(idx_overlaps) > 1:
            replacement_required = True
            idx_overlaps.remove(i)  # don't move yourself
            # dist is mostly negative, because of overlap
            replace_size = min_dist_between_objects - dist[idx_overlaps]
            xy = radial_replacement_from_reference_dots(
                xy=xy,
                ref_pos_id=i,
                neighbour_ids=idx_overlaps,
                replacement_size=replace_size)
    return replacement_required


class BrownianMotion(object):

    def __init__(self, start_pos, delta=2, search_area_radius=None, bounce_boarder=True):
        """performs brownian motions (search walk) optionally in a circular area

        Parameters
        ----------

        bounce_boarder: bool
            if true, random walk bounces back at circle boarder, otherwise walk
            will be continued in the center of the area.

        Notes
        -----
        see Brownian motion https://en.wikipedia.org/wiki/Brownian_motion
        """
        if search_area_radius is not None and \
                np.hypot(start_pos[0], start_pos[1]) > search_area_radius:
            raise ValueError("start_pos is outside max_radius")

        self.area_radius = search_area_radius
        self.scale = delta ** 2
        self._history = [np.array(start_pos)]
        self.bounce = bounce_boarder

    @property
    def current(self):
        return self._history[-1]

    def step_back(self):
        """redo last step """
        if len(self._history) > 1:
            return self._history.pop()

    def center_last_step(self):
        last = self.step_back()
        if last is not None:
            new = self._history[-1] - last  # type: ignore
            self._history.append(new)
        return self.current

    def next(self, dt=1):
        while True:
            new = rng.generator.normal(loc=0, scale=self.scale * dt,
                                       size=2) + self.current
            if self.area_radius is None or \
                    np.hypot(new[0], new[1]) <= self.area_radius:
                self._history.append(new)
                return self.current

            elif not self.bounce:
                new = new - self.current  # do this step from center
                self._history.append(new)
                return self.current
