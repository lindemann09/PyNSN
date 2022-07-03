import random
import numpy as np

from .._lib import  geometry


def _jitter_identical_positions(xy, jitter_size=0.1):
    """jitters points with identical position"""

    for idx, ref_object in enumerate(xy):
        identical = np.where(np.all(np.equal(xy, ref_object), axis=1))[0]  # find identical positions
        if len(identical) > 1:
            for x in identical:  # jitter all identical positions
                if x != idx:
                    xy[x, :] = xy[x, :] - geometry.polar2cartesian(
                        [[jitter_size, random.random() * 2 * np.pi]])[0]

    return xy

def radial_replacement_from_reference_dots(xy, ref_pos_id,
                                            neighbour_ids, replacement_size):
    """remove neighbouring position radially from reference position
    helper function, typically used for realign
    """

    # check if there is an identical position and jitter to avoid fully overlapping positions
    if np.sum(np.all(xy[neighbour_ids,] == xy[ref_pos_id, :],
                   axis=1)) > 0:
        xy = _jitter_identical_positions(xy)

    # relative polar positions to reference_dot
    tmp_polar = geometry.cartesian2polar(xy[neighbour_ids, :] - xy[ref_pos_id, :])
    tmp_polar[:, 0] = 0.000000001 + replacement_size # determine movement size
    tmp_xy = geometry.polar2cartesian(tmp_polar)
    xy[neighbour_ids, :] = np.array([xy[neighbour_ids, 0] + tmp_xy[:, 0],
                                           xy[neighbour_ids, 1] + tmp_xy[:, 1]]).T
    return xy

def remove_overlap_from_inner_to_outer(xy, min_dist_between, distance_matrix_function):
    """returns xy and boolean, if replacements were required"""
    assert callable(distance_matrix_function)

    replacement_required = False
    # from inner to outer remove overlaps
    for i in np.argsort(geometry.cartesian2polar(xy, radii_only=True)):
        dist_mtx = distance_matrix_function(between_positions=False)
        dist = dist_mtx[i,:]
        idx_overlaps = np.where(dist < min_dist_between)[0].tolist()  # overlapping dot ids
        if len(idx_overlaps) > 1:
            replacement_required = True
            idx_overlaps.remove(i)  # don't move yourself
            replace_size =min_dist_between - dist[idx_overlaps]  # dist is mostly negative, because of overlap
            xy = radial_replacement_from_reference_dots(
                                         xy=xy,
                                         ref_pos_id=i,
                                         neighbour_ids=idx_overlaps,
                                         replacement_size=replace_size)
    return replacement_required
