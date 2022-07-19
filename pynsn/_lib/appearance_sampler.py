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


class AppearanceSampler(object):

    def __init__(self):

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

            return [shapes.Dot(xy=(0, 0), diameter=dia, attribute=attr) \
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

            return [shapes.Rectangle(xy=(0, 0), size=(w,h), attribute=attr)\
                    for w, h, attr in zip(width, height, attributes)]

    def as_dict(self):
        rtn = {}
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
