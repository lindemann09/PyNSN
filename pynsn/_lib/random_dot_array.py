__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'


from . import _misc as _misc
from ._item_attributes import ItemAttributes as _ItemAttributes
from ._dot_array import DotArray as _DotArray


class Specs(object):

    def __init__(self,
                 target_area_radius,
                 item_diameter_mean,
                 item_diameter_range=None,
                 item_diameter_std=None,
                 item_colour=None,  # todo feature
                 minimum_gap=2,
                 min_distance_area_boarder=0):
        # minim gap

        """Specification of a Random Dot Array

        Parameters:
        -----------
        stimulus_area_radius : int
            the radius of the stimulus area
        n_dots : int
            number of moving dots

        automatic logging log only the create process. If colours a changes later they are not log.
        Use manual logging in this case.

        """

        if item_diameter_std <= 0:
            item_diameter_std = None
        elif item_diameter_range is not None and \
                (item_diameter_mean <= item_diameter_range[0] or
                 item_diameter_mean >= item_diameter_range[1] or
                 item_diameter_range[0] >= item_diameter_range[1]):
            raise RuntimeError("item_diameter_mean has to be inside the defined item_diameter_range")

        self.minimum_gap = minimum_gap
        self.target_array_radius = target_area_radius
        self.item_diameter_range = item_diameter_range
        self.item_diameter_mean = item_diameter_mean
        self.item_diameter_std = item_diameter_std
        self.min_distance_area_boarder = min_distance_area_boarder
        self.item_attributes = _ItemAttributes(colour=item_colour)

    def as_dict(self):
        return {"target_array_radius": self.target_array_radius,
                "dot_diameter_mean": self.item_diameter_mean,
                "dot_diameter_range": self.item_diameter_range,
                "dot_diameter_std": self.item_diameter_std,
                "dot_colour": str(self.item_attributes.colour),  ##todo feature
                "minimum_gap": self.minimum_gap,
                "min_distance_area_boarder": self.min_distance_area_boarder}


def create(n_dots, specs, occupied_space=None,
           logger=None):
    """occupied_space is a dot array (used for multicolour dot array (join after)

    returns None if not possible
    """

    assert isinstance(specs, Specs)

    rtn = _DotArray(target_array_radius=specs.target_array_radius,
                   minimum_gap=specs.minimum_gap)

    # random diameter from beta distribution with exact mean and str
    diameters = _misc.random_beta(size=n_dots,
                                 number_range=specs.item_diameter_range,
                                 mean=specs.item_diameter_mean,
                                 std=specs.item_diameter_std)

    for dia in diameters:
        try:
            xy = rtn.random_free_dot_position(dot_diameter=dia,
                                              occupied_space=occupied_space,
                                              min_distance_area_boarder=
                                            specs.min_distance_area_boarder)
        except:
            return None
        rtn.append(xy=xy, item_diameters=dia,
                     attributes=specs.item_attributes)

    if logger is not None:
        from .logging import LogFile # to avoid circular import
        if not isinstance(logger, LogFile): #
            raise RuntimeError("logger has to be None or a GeneratorLogger")
        logger.log(rtn)

    return rtn