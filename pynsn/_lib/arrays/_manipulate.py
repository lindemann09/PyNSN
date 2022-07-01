import random
from copy import deepcopy
import numpy as np
from scipy import spatial

from .. import  geometry
from .. import shapes


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
        dist_mtx = distance_matrix_function(between_positions=False,
                                         overlap_is_zero=False)
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


def get_random_free_position(the_object,
                             target_area_radius,
                             allow_overlapping,
                             distances_function,
                             min_dist_between,
                             min_dist_area_boarder,
                             occupied_space_distances_function,
                             convex_hull_xy):
    """returns a available random xy position

    raise exception if not found
    occupied space: see generator generate
    """

    N_ATTEMPTS = 3000
    assert callable(distances_function)
    if isinstance(the_object, shapes.Dot):
        object_size = the_object.diameter / 2.0
    elif isinstance(the_object, shapes.Rectangle):
        object_size = max(the_object.size)
    else:
        raise NotImplementedError()
    inside_convex_hull = convex_hull_xy is not None
    if inside_convex_hull:
        delaunay = spatial.Delaunay(convex_hull_xy)
    else:
        delaunay = None

    area_rad = target_area_radius - min_dist_area_boarder - object_size
    proposal_obj = deepcopy(the_object)
    cnt = 0
    while True:
        cnt += 1
        ##  polar method seems to produce central clustering
        #  proposal_polar =  np.array([random.random(), random.random()]) *
        #                      (target_radius, TWO_PI)
        proposal_obj.xy = np.array([random.random(), random.random()]) \
                           * 2 * area_rad - area_rad

        # is outside area
        if isinstance(the_object, shapes.Dot):
            bad_position = area_rad <= proposal_obj.polar_radius
        else:
            # Rect: check if one edge is outside
            bad_position = False
            for e in proposal_obj.edges():
                if e.polar_radius >= area_rad:
                    bad_position = True
                    break

        if not bad_position and inside_convex_hull:
            bad_position = delaunay.find_simplex(np.array(proposal_obj.xy)) < 0

        if not bad_position and not allow_overlapping:
            # find bad_positions
            dist = distances_function(proposal_obj)
            if callable(occupied_space_distances_function):
                dist = np.append(dist, occupied_space_distances_function(proposal_obj))
            idx = np.where(dist < min_dist_between)[0]  # overlapping dot ids
            bad_position = len(idx) > 0

        if not bad_position:
            return proposal_obj.xy
        elif cnt > N_ATTEMPTS:
            raise StopIteration(u"Can't find a free position")
