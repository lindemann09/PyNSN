__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import copy
from random import shuffle

from .. import shapes
from .. import arrays
from ._size_distribution import SizeDistribution

def create(reference_array,
           size_distribution,
           n_objects,
           attributes = None,
           allow_overlapping = False,
           occupied_space = None):
    """occupied_space is a dot array (used for multicolour dot array (join after)

    attribute is an array, arrays are assigned randomly.

    """
    if not isinstance(reference_array, arrays.GenericObjectArray):
        raise TypeError("Reference array has to be of type DotArray, RectangleArray or GenericObjectArray, but not {}".format(
                        type(reference_array).__name__))
    if not isinstance(size_distribution, SizeDistribution):
        raise RuntimeError("Size distribution has to be of type SizeDistribution, but not {}".format(
                        type(size_distribution).__name__))


    if size_distribution.diameter is not None:
        # DotArray
        rtn = arrays.DotArray(target_area_radius=reference_array.target_area_radius,
                       min_dist_between=reference_array.min_dist_between,
                       min_dist_area_boarder=reference_array.min_dist_area_boarder)

        for dot in size_distribution.sample(n=n_objects):
            try:
                dot = rtn.get_random_free_position(ref_object=dot,
                                                   occupied_space=occupied_space,
                                                   allow_overlapping=allow_overlapping)
            except StopIteration as e:
                raise StopIteration("Can't find a solution for {} items in this array".format(n_objects))
            rtn.add([dot])

    else:
        # RectArray
        rtn = arrays.RectangleArray(target_area_radius=reference_array.target_area_radius,
                             min_dist_between=reference_array.min_dist_between,
                             min_dist_area_boarder=reference_array.min_dist_area_boarder)

        for rect in size_distribution.sample(n=n_objects):
            try:
                rect = rtn.get_random_free_position(ref_object=rect,
                                                    occupied_space=occupied_space,
                                                    allow_overlapping=allow_overlapping)
            except StopIteration:
                raise StopIteration("Can't find a solution for {} items in this array".format(n_objects))

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