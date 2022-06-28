from collections import OrderedDict
from ..misc import dict_to_text

from ..distributions import PyNSNDistribution

class SizeDistribution(object):

    def __init__(self, diameter=None, width=None, height=None):

         is_rect = width is not None or height is not None
         is_dot = diameter is not None
         if is_dot and is_rect:
             raise TypeError("Please define either diameter or width and height, not both.")
         elif is_rect and (width is None or height is None):
             raise TypeError("Please define width and height for rectangles.")
         elif not is_dot and not is_rect:
             raise TypeError("No size distribution define. Please define either diameter or width and height.")

         if diameter is not None and not isinstance(diameter, PyNSNDistribution):
            raise TypeError("diameter has to be a PyNSNDistribution")
         if width is not None and not isinstance(width, PyNSNDistribution):
            raise TypeError("width has to be a PyNSNDistribution")
         if height is not None and not isinstance(height, PyNSNDistribution):
            raise TypeError("height has to be a PyNSNDistribution")


         self.diameter = diameter
         self.width = width
         self.height = height

    def as_dict(self):
        try:
            d = self.diameter.as_dict()
        except:
            d = None
        try:
            w = self.width.as_dict()
        except:
            w = None
        try:
            h = self.height.as_dict()
        except:
            h = None
        return OrderedDict({
            "diameter": d, "width": w, "height": h})

    def __str__(self):
        return dict_to_text(self.as_dict(), col_a=12, col_b=7)


