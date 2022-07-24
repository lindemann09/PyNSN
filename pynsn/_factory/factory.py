__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import copy

from .._arrays.dot_array import DotArray
from .._arrays.parameter import ArrayParameter
from .._arrays.rect_array import RectangleArray
from .._shapes.dot import Dot
from .._shapes.rectangle import Rectangle
from .._lib.misc import dict_to_text
from .distributions import PyNSNDistribution, _round_samples, Levels
from .._lib.exception import NoSolutionError
from .._lib.lib_typing import Union

# FIXME typing required
class _Constant(object):

    def __init__(self, value):
        """Helper class to "sample" constance.

        Looks like a PyNSNDistribution, but sample returns just the constant

        Parameter:
        ----------
        constant : numeric
        """

        self.value = value

    def sample(self, n, round_to_decimals=None):
        return _round_samples([self.value] * n, round_to_decimals)

    def as_dict(self):
        return {"distribution": "Constant",
                "value": self.value}


def _make_distr(value):
    """helper
    returns a distribution or None, if None
    """

    if value is None:
        return None
    elif isinstance(value, PyNSNDistribution):
        return value
    elif isinstance(value, (list, tuple)):
        return Levels(levels=copy(value))
    else:
        return _Constant(value)


# TODO typing


class NSNFactory(ArrayParameter):

    def __init__(self, target_area_radius,
                 min_dist_between=None,
                 min_dist_area_boarder=None):
        """

        Parameters
        ----------
        target_area_radius
        min_dist_between
        min_dist_area_boarder
        """

        ArrayParameter.__init__(self,
                                target_area_radius=target_area_radius,
                                min_dist_between=min_dist_between,
                                min_dist_area_boarder=min_dist_area_boarder)
        self._distr_diameter = None
        self._distr_width = None
        self._distr_height = None
        self._distr_proportion = None
        self._distr_attributes = None

    @property
    def distr_diameter(self):
        return self._distr_diameter

    @property
    def distr_width(self):
        return self._distr_width

    @property
    def distr_height(self):
        return self._distr_height

    @property
    def distr_proportion(self):
        return self._distr_proportion

    @property
    def distr_attributes(self):
        return self._distr_attributes

    def set_appearance_dot(self, diameter, attributes=None):
        self._distr_width = None
        self._distr_height = None
        self._distr_proportion = None
        self._distr_diameter = _make_distr(diameter)
        self._distr_attributes = _make_distr(attributes)

    def set_appearance_rectangle(self, width=None, height=None,
                                 proportion=None, attributes=None):

        n_rect_parameter = sum([width is not None, height is not None,
                                proportion is not None])
        if n_rect_parameter == 1:
            raise TypeError("Define rectangle width and height or, alternatively, rectangle proportion together with "
                            "either width or height.")
        self._distr_diameter = None
        self._distr_width = _make_distr(width)
        self._distr_height = _make_distr(height)
        self._distr_proportion = _make_distr(proportion)
        self._distr_attributes = _make_distr(attributes)

    def is_appearance_set(self):
        return self._distr_diameter is not None or self._distr_width is not None or \
               self._distr_height is not None

    def sample(self, n, round_to_decimals=None):
        """return list objects (Dot or Rect) with random size
        all positions = (0,0)
        """
        if self._distr_attributes is not None:
            attributes = self._distr_attributes.sample(n)
        else:
            attributes = [None] * n
        if self._distr_diameter is not None:
            diameter = self._distr_diameter.sample(n)

            return [Dot(xy=(0, 0), diameter=dia, attribute=attr)
                    for dia, attr in zip(diameter, attributes)]
        else:
            # Rect
            try:
                width = self._distr_width.sample(n)
            except AttributeError:
                width = None
            try:
                height = self._distr_height.sample(n)
            except AttributeError:
                height = None
            try:
                proportion = self._distr_proportion.sample(n)
            except AttributeError:
                proportion = None

            if height is None:
                height = width * proportion
            elif width is None:
                width = height / proportion

            if round_to_decimals is not None:
                width = _round_samples(width, round_to_decimals=round_to_decimals)
                height = _round_samples(width, round_to_decimals=round_to_decimals)

            return [Rectangle(xy=(0, 0), size=(w, h), attribute=attr)
                    for w, h, attr in zip(width, height, attributes)]

    def as_dict(self):
        rtn = ArrayParameter.as_dict(self)
        try:
            rtn.update({"diameter": self._distr_diameter.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"width": self._distr_width.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"height": self._distr_height.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"proportion": self._distr_proportion.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"attributes": self._distr_attributes.as_dict()})
        except AttributeError:
            pass
        return rtn

    def create_random_array(self, n_objects: int,
                            allow_overlapping: bool =False,
                            occupied_space: Union[None, DotArray, RectangleArray] = None) -> Union[DotArray, RectangleArray]:
        """
        occupied_space is a dot array (used for multicolour dot array (join after)

        attribute is an array, _arrays are assigned randomly.


        Parameters
        ----------
        n_objects
        allow_overlapping
        occupied_space

        Returns
        -------
        rtn : object array
        """
        if not self.is_appearance_set:
            raise RuntimeError("No appearance defined. Please use 'set_dot', 'set_rect_or' "
                               "'set_appearance'")
        if self._distr_diameter is not None:
            # DotArray
            rtn = DotArray(target_area_radius=self.target_area_radius,
                           min_dist_between=self.min_dist_between,
                           min_dist_area_boarder=self.min_dist_area_boarder)

            for dot in self.sample(n=n_objects):
                try:
                    dot = rtn.get_free_position(ref_object=dot, in_neighborhood=False,
                                                occupied_space=occupied_space,
                                                allow_overlapping=allow_overlapping)
                except NoSolutionError:
                    raise NoSolutionError("Can't find a solution for {} items in this array".format(n_objects))
                rtn.add([dot])

        else:
            # RectArray
            rtn = RectangleArray(target_area_radius=self.target_area_radius,
                                 min_dist_between=self.min_dist_between,
                                 min_dist_area_boarder=self.min_dist_area_boarder)

            for rect in self.sample(n=n_objects):
                try:
                    rect = rtn.get_free_position(ref_object=rect, in_neighborhood=False,
                                                 occupied_space=occupied_space,
                                                 allow_overlapping=allow_overlapping)
                except NoSolutionError:
                    raise NoSolutionError("Can't find a solution for {} ".format(n_objects) +
                                          "items in this array.")

                rtn.add([rect])
        return rtn

    def create_incremental_random_array(self, n_objects,
                                        allow_overlapping=False):
        """

        Parameters
        ----------
        n_objects
        allow_overlapping

        Returns
        -------
        rtn : iterator of object _arrays
        """
        previous = None
        for n in range(n_objects):
            current = self.create_random_array(n_objects=1,
                                               allow_overlapping=allow_overlapping,
                                               occupied_space=previous)
            if previous is not None:
                current.join(previous)
            previous = current
            yield current

    def __str__(self):
        return dict_to_text(self.as_dict(), col_a=12, col_b=7)
