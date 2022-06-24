__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import copy as _copy
from random import shuffle as _shuffle

from ._lib import distributions as distr
from ._lib import shape as _shape
from ._lib.arrays import DotArray as _DotArray
from ._lib.arrays import RectangleArray as _RectangleArray


class _Specs(object):

    def __init__(self,  target_area_radius,
                 minimum_gap,
                 min_distance_area_boarder):
        self.minimum_gap = minimum_gap
        self.target_array_radius = target_area_radius
        self.min_distance_area_boarder = min_distance_area_boarder

    def as_dict(self):
        return {"target_array_radius": self.target_array_radius,
                "minimum_gap": self.minimum_gap,
                "min_distance_area_boarder": self.min_distance_area_boarder}

    def copy(self):
        """returns a deepcopy of the specs"""
        return _copy.deepcopy(self)


class DotArraySpecs(_Specs):

    def __init__(self,
                 target_area_radius,
                 diameter_distribution,
                 minimum_gap = 2,
                 min_distance_area_boarder = 1):
        """

        Parameters
        ----------
        target_area_radius
        diameter_distribution
        minimum_gap
        min_distance_area_boarder
        """
        super().__init__(target_area_radius=target_area_radius,
                         minimum_gap=minimum_gap,
                         min_distance_area_boarder=min_distance_area_boarder)
        if not isinstance(diameter_distribution, distr.PyNSNDistribution):
            raise TypeError("diameter_distribution has to be a PyNSNDistribution")
        self.diameter_distr = diameter_distribution

    def as_dict(self):
        rtn = super().as_dict()
        rtn.update({"diameter_distr": self.diameter_distr.as_dict()})
        return rtn


class RectangleArraySpecs(_Specs):

    def __init__(self,
                 target_area_radius,
                 width_distribution,
                 height_distribution,
                 minimum_gap = 2,
                 min_distance_area_boarder = 1):
        """

        Parameters
        ----------
        target_area_radius
        width_distribution
        height_distribution
        minimum_gap
        """
        super().__init__(target_area_radius=target_area_radius,
                         minimum_gap=minimum_gap,
                         min_distance_area_boarder=min_distance_area_boarder)
        if not isinstance(width_distribution, distr.PyNSNDistribution):
            raise TypeError("width_distribution has to be a PyNSNDistribution")
        if not isinstance(height_distribution, distr.PyNSNDistribution):
            raise TypeError("height_distribution has to be a PyNSNDistribution")
        self.width_distr = width_distribution
        self.height_distr = height_distribution

    def as_dict(self):
        rtn = super().as_dict()
        rtn.update({"width_distr": self.width_distr.as_dict(),
                    "height_distr": self.height_distr.as_dict()})
        return rtn


def random_array(specs,
                 n_objects,
                 attributes = None,
                 allow_overlapping = False,
                 occupied_space = None):
    """occupied_space is a dot array (used for multicolour dot array (join after)

    attribute is an array, arrays are assigned randomly.

    """

    if isinstance(specs, DotArraySpecs):
        # DotArray
        rtn = _DotArray(target_array_radius=specs.target_array_radius,
                        minimum_gap=specs.minimum_gap)

        for dia in specs.diameter_distr.sample(n=n_objects):
            try:
                xy = rtn.random_free_position(dot_diameter=dia,
                          occupied_space=occupied_space,
                          allow_overlapping=allow_overlapping,
                          min_distance_area_boarder=specs.min_distance_area_boarder)
            except StopIteration as e:
                raise StopIteration("Can't find a solution for {} items this dot".format(n_objects))
            rtn.add([_shape.Dot(xy=xy, diameter=dia)])

    elif isinstance(specs, RectangleArraySpecs):
        # RectArray
        rtn = _RectangleArray(target_array_radius=specs.target_array_radius,
                        minimum_gap=specs.minimum_gap)

        sizes = zip(specs.width_distr.sample(n=n_objects),
                    specs.height_distr.sample(n=n_objects))

        for s in sizes:
            try:
                xy = rtn.random_free_position(rectangle_size=s,
                          occupied_space=occupied_space,
                          allow_overlapping=allow_overlapping,
                          min_distance_area_boarder=specs.min_distance_area_boarder)
            except StopIteration as e:
                raise StopIteration("Can't find a solution for {} items this dot".format(n_objects))

            rtn.add([_shape.Rectangle(xy=xy, size=s)])

    else:
        raise RuntimeError("specs has to be of type DotArraySpecs or , but not {}".format(
                        type(specs).__name__))

    # attribute assignment
    if isinstance(attributes, (tuple, list)):
        att = []
        while len(att) < n_objects:
            tmp = _copy.copy(attributes)
            _shuffle(tmp)
            att.extend(tmp)
        _shuffle(att)
        rtn.set_attributes(att[:n_objects])
    else:
        rtn.set_attributes(attributes)

    return rtn