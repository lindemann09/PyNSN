__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import copy
from random import shuffle

from .._lib import DotArray, RectangleArray, _check_array_parameter, NoSolutionError
from ._size_distribution import SizeDistribution


def create(array_parameter,
           size_distribution,
           n_objects,
           attributes = None,
           allow_overlapping = False,
           occupied_space = None):
    """occupied_space is a dot array (used for multicolour dot array (join after)

    attribute is an array, arrays are assigned randomly.

    """
    _check_array_parameter(array_parameter)
    if not isinstance(size_distribution, SizeDistribution):
        raise RuntimeError("Size distribution has to be of type SizeDistribution, but not {}".format(
                        type(size_distribution).__name__))

    if size_distribution.diameter is not None:
        # DotArray
        rtn = DotArray(target_area_radius=array_parameter.target_area_radius,
                       min_dist_between=array_parameter.min_dist_between,
                       min_dist_area_boarder=array_parameter.min_dist_area_boarder)

        for dot in size_distribution.sample(n=n_objects):
            try:
                dot = rtn.get_random_free_position(ref_object=dot,
                                                   occupied_space=occupied_space,
                                                   allow_overlapping=allow_overlapping)
            except NoSolutionError as e:
                raise NoSolutionError("Can't find a solution for {} items in this array".format(n_objects))
            rtn.add([dot])

    else:
        # RectArray
        rtn = RectangleArray(target_area_radius=array_parameter.target_area_radius,
                             min_dist_between=array_parameter.min_dist_between,
                             min_dist_area_boarder=array_parameter.min_dist_area_boarder)

        for rect in size_distribution.sample(n=n_objects):
            try:
                rect = rtn.get_random_free_position(ref_object=rect,
                                                    occupied_space=occupied_space,
                                                    allow_overlapping=allow_overlapping)
            except NoSolutionError:
                raise NoSolutionError("Can't find a solution for {} items in this array".format(n_objects))

            rtn.add([rect])

    # attribute assignment
    if isinstance(attributes, (tuple, list)):
        att = []
        while len(att) < n_objects:
            tmp = copy.copy(attributes)
            shuffle(tmp)
            att.extend(tmp)
        shuffle(att)
        rtn.set_attributes(att[:n_objects])
    else:
        rtn.set_attributes(attributes)

    return rtn


def create_incremental(reference_array,
                       size_distribution,
                       n_objects,
                       attributes = None,
                       allow_overlapping = False):

    previous = None
    for n in range(n_objects):
        current = create(array_parameter=reference_array,
                         size_distribution=size_distribution,
                         n_objects=1,
                         attributes=attributes,
                         allow_overlapping=allow_overlapping,
                         occupied_space=previous)
        if previous is not None:
            current.join(previous)
        previous = current
        yield current

