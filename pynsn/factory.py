__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import copy as _copy
from ._nsn.dot_array import DotArray as _DotArray
from ._nsn.rect_array import RectangleArray as _RectangleArray
from ._nsn import shape as _shape
from ._lib import distributions as distr


class _Specs(object):

    def __init__(self, target_area_radius, minimum_gap,
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
                 minimum_gap=2,
                 min_distance_area_boarder=1):
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
        if not isinstance(diameter_distribution, distr._PyNSNDistribution):
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
                 minimum_gap=2,
                 min_distance_area_boarder=1):
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
        if not isinstance(width_distribution, distr._PyNSNDistribution):
                raise TypeError("width_distribution has to be a PyNSNDistribution")
        if not isinstance(height_distribution, distr._PyNSNDistribution):
                raise TypeError("height_distribution has to be a PyNSNDistribution")
        self.width_distr = width_distribution
        self.height_distr = height_distribution

    def as_dict(self):
        rtn = super().as_dict()
        rtn.update({"width_distr": self.width_distr.as_dict(),
                    "height_distr": self.height_distr.as_dict()})
        return rtn


def random_array(specs, n_dots, attribute=None, occupied_space=None):
    """occupied_space is a dot array (used for multicolour dot array (join after)

    returns None if not possible
    """

    if isinstance(specs, DotArraySpecs):
        # DotArray
        rtn = _DotArray(target_array_radius=specs.target_array_radius,
                        minimum_gap=specs.minimum_gap)

        for dia in specs.diameter_distr.sample(n=n_dots):
            try:
                xy = rtn.random_free_position(dot_diameter=dia,
                          occupied_space=occupied_space,
                          min_distance_area_boarder=specs.min_distance_area_boarder)
            except:
                return None
            rtn.add([_shape.Dot(xy=xy, diameter=dia, attribute=attribute)])

    elif isinstance(specs, RectangleArraySpecs):
        # RectArray
        rtn = _RectangleArray(target_array_radius=specs.target_array_radius,
                        minimum_gap=specs.minimum_gap)

        sizes = zip(specs.width_distr.sample(n=n_dots),
                    specs.height_distr.sample(n=n_dots))
        for s in sizes:
            try:
                xy = rtn.random_free_position(rectangle_size=s,
                      occupied_space=occupied_space,
                      min_distance_area_boarder=specs.min_distance_area_boarder)
            except:
                return None
            rtn.add([_shape.Rectangle(xy=xy, size=s, attribute=attribute)])

    else:
        raise TypeError("specs has to be of type DotArraySpecs or , but not {}".format(
                        type(specs).__name__))
    return rtn