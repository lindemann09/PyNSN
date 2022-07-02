from .._lib.misc import dict_to_text
from .. import shapes
from ..distributions import PyNSNDistribution, Constant, _round_samples


class SizeDistribution(object):

    def __init__(self, diameter=None, width=None,
                 height=None, rectangle_proportion=None):

        rect_parameter_set = [width is not None,
                              height is not None,
                              rectangle_proportion is not None]

        is_rect = sum(rect_parameter_set)>0
        is_dot = diameter is not None
        if is_dot and is_rect:
         raise TypeError("Please define either dot or rectangle parameter, not both.")
        elif not is_dot and not is_rect:
         raise TypeError("No size distribution define. Please define either diameter or width and height.")
        elif is_rect and sum(rect_parameter_set) != 2:
            raise TypeError("Define rectangle width and height or, alternatively, rectangle proportion together with either width or height.")

        if isinstance(diameter, (int, float)):
            diameter = Constant(diameter)
        if isinstance(width, (int, float)):
            width = Constant(width)
        if isinstance(height, (int, float)):
            height = Constant(height)
        if isinstance(rectangle_proportion, (int, float)):
            rectangle_proportion = Constant(rectangle_proportion)

        err = "{} has to be a PyNSNDistribution or single number (constant)"
        if diameter is not None and not isinstance(diameter, PyNSNDistribution):
            raise TypeError(err.format("dot_diameter"))
        if width is not None and not isinstance(width, PyNSNDistribution):
            raise TypeError(err.format("rectangle_width"))
        if height is not None and not isinstance(height, PyNSNDistribution):
            raise TypeError(err.format("rectangle_height"))
        if rectangle_proportion is not None and \
                not isinstance(rectangle_proportion, PyNSNDistribution):
            raise TypeError(err.format("rectangle_proportion"))

        self.diameter = diameter
        self.width = width
        self.height = height
        self.proportion = rectangle_proportion

    def sample(self, n, round_to_decimals=None):
        """return list objects (Dot or Rect) with random size
        all positions = (0,0)
        """
        if self.diameter is not None:
            diameter = self.diameter.sample(n)
            return [shapes.Dot(xy=(0, 0), diameter=dia) for dia in diameter]
        else:
            #Rect
            try:
                width = self.width.sample(n)
            except AttributeError:
                width = None
            try:
                height = self.height.sample(n)
            except AttributeError:
                height = None
            try:
                proportion = self.proportion.sample(n)
            except AttributeError:
                proportion = None

            if height is None:
                height = width * proportion
            elif width is None:
                width = height / proportion

            if round_to_decimals is not None:
                width = _round_samples(width, round_to_decimals=round_to_decimals)
                height = _round_samples(width, round_to_decimals=round_to_decimals)

            return [shapes.Rectangle(xy=(0, 0), size=s) for s in zip(width, height)]


    def as_dict(self):
        rtn = {}
        try:
            rtn.update({"diameter": self.diameter.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"width": self.width.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"height": self.height.as_dict()})
        except AttributeError:
            pass
        try:
            rtn.update({"proportion": self.proportion.as_dict()})
        except AttributeError:
            pass
        return rtn


    def __str__(self):
        return dict_to_text(self.as_dict(), col_a=12, col_b=7)


