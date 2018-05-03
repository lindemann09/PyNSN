from __future__ import division

from math import sqrt, pi
try:
    from math import log2
except:
    from math import log
    log2 = lambda x: log(x, 2)

class CardinalFeatures(object):

    def __init__(self, Numerosity=1, Size=1, Spacing=1):
        """
        :param Numerosity:
        :param Size:
        :param Spacing:
        """

        self.Numerosity = Numerosity
        self.Size = Size
        self.Spacing = Spacing

    # -- cardinal features

    @property
    def Numerosity(self):
        return self._N

    @Numerosity.setter
    def Numerosity(self, x):
        self._N = x
        self._sqrtN = sqrt(x)

    @property
    def Size(self):
        return self._size

    @Size.setter
    def Size(self, x):
        self._size = float(x)
        self._sqrtSize = sqrt(x)

    @property
    def Spacing(self):
        return self._spacing

    @Spacing.setter
    def Spacing(self, x):
        self._spacing = float(x)
        self._sqrtSpace = sqrt(x)


    # -- stimulus features

    @property
    def total_surface_area(self):
        """TSA = sqrt(Sz*n) """
        return self._sqrtSize * self._sqrtN

    @total_surface_area.setter
    def total_surface_area(self, x):
        self.Size = (x / self._sqrtN) ** 2

    @property
    def item_surface_area(self):
        """mean item surface array. ISA = sqrt(Size/N)"""
        return self._sqrtSize / self._sqrtN

    @item_surface_area.setter
    def item_surface_area(self, x):
        """mean item surface array. ISA = sqrt(Size/N)"""
        self.Size = (x * self._sqrtN) ** 2

    @property
    def field_area(self):
        """FA = sqrt(Spacing * N)"""
        return self._sqrtSpace * self._sqrtN

    @field_area.setter
    def field_area(self, x):
        self.Spacing = (x / self._sqrtN) ** 2

    @property
    def sparsity(self):
        return self._sqrtSpace / self._sqrtN

    @sparsity.setter
    def sparsity(self, x):
        self.Spacing = (x * self._sqrtN) ** 2

    # ---- further stimulus features (only setter)
    @property
    def density(self):
        return 1/float(self.sparsity)

    @property
    def coverage(self):
        """often referred to as density: total surface area per field area

        Can not be modified directly, because it reflects relation of two cardinal feature.
        Please adjust Size and/or Spacing.
        """
        return self._sqrtSize / self._sqrtSpace

    @property
    def apparent_closeness(self):
        """

        Can not be modified directly, because it reflects relation of two cardinal feature.
        Please adjust Size and/or Spacing.
        """

        return self._sqrtSize * self._sqrtSpace

class CardinalFeaturesDotArray(CardinalFeatures):

    @property
    def item_diameter(self):

        return sqrt(self.item_surface_area/pi)  * 2

    @item_diameter.setter
    def item_diameter(self, x):
        self.item_surface_area = pi * x**2 / 4

    @property
    def item_perimeter(self):
        return 2 * sqrt(pi * self.item_surface_area)

    @item_perimeter.setter
    def item_perimeter(self, x):
        self.item_surface_area = x**2 / (4*pi)


