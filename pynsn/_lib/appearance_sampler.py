from copy import copy

from .._lib import shapes
from ..distributions import PyNSNDistribution, _round_samples, Levels


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


def _type_check_distr(value, variable_name):
    """helper
    returns a distribution
    """

    #err = "{} has to be a PyNSNDistribution, a single number (constant) or None"
    #if value is not None and not isinstance(value, (int, float, PyNSNDistribution)):
    #    raise TypeError(err.format(variable_name))

    if value is None:
        return None
    elif isinstance(value, PyNSNDistribution):
        return value
    elif isinstance(value, (list, tuple)):
        return Levels(levels=copy(value))
    else:
        return _Constant(value)


class AppearanceSampler(object):

    def __init__(self):

        self._diameter = None
        self._width = None
        self._height = None
        self._proportion = None
        self._attributes = None

    def set_appearance_dot(self, diameter, attributes=None):
        self.set_appearance_rectangle(width=None, height=None, proportion=None)
        self._diameter = _type_check_distr(diameter, "diameter")
        self._attributes = _type_check_distr(attributes, "attributes")

    def set_appearance_rectangle(self, width=None, height=None,
                                 proportion=None, attributes=None):
        n_rect_parameter = sum([width is not None, height is not None,
                                proportion is not None])
        if n_rect_parameter == 1:
            raise TypeError("Define rectangle width and height or, alternatively, rectangle proportion together with "
                            "either width or height.")
        self._width = _type_check_distr(width, "width")
        self._height = _type_check_distr(height, "height")
        self._proportion = _type_check_distr(proportion, "proportion")
        self._attributes = _type_check_distr(attributes, "attributes")

    def is_appearance_set(self):
        return self._diameter is not None or self._width is not None or \
               self._height is not None

    def sample(self, n, round_to_decimals=None):
        """return list objects (Dot or Rect) with random size
        all positions = (0,0)
        """
        if self._attributes is not None:
            attributes = self._attributes.sample(n)
        else:
            attributes = [None] * n
        if self._diameter is not None:
            diameter = self._diameter.sample(n)

            return [shapes.Dot(xy=(0, 0), diameter=dia, attribute=attr) \
                    for dia, attr in zip(diameter, attributes)]
        else:
            # Rect
            try:
                width = self._width.sample(n)
            except AttributeError:
                width = None
            try:
                height = self._height.sample(n)
            except AttributeError:
                height = None
            try:
                proportion = self._proportion.sample(n)
            except AttributeError:
                proportion = None

            if height is None:
                height = width * proportion
            elif width is None:
                width = height / proportion

            if round_to_decimals is not None:
                width = _round_samples(width, round_to_decimals=round_to_decimals)
                height = _round_samples(width, round_to_decimals=round_to_decimals)

            return [shapes.Rectangle(xy=(0, 0), size=(w,h), attribute=attr)\
                    for w, h, attr in zip(width, height, attributes)]

    def as_dict(self):
        rtn = {}
        try:
            rtn.update({"diameter": self._diameter.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"width": self._width.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"height": self._height.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"proportion": self._proportion.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"attributes": self._attributes.as_dict()})
        except AttributeError:
            pass
        return rtn
