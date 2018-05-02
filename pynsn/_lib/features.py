from math import log2



class ArrayFeatures(object):

    def __init__(self, Numerosity,
                 total_surface_area=None,
                 item_surface_area=None,
                 field_area=None,
                 sparcity=None,
                 Size=None,
                 Spacing=None):
        """
        [N, [TSA, ISA], [FA, SPAR]]

        Please define two features for two diffienct categories (A, B C).


        A: TSA (total surface area), ISA (mean item surface area), Size
        B: FA (Field area), SPAR (Sparcity), Spacing
        """

        self.N = Numerosity

        self.TSA = total_surface_area
        self.ISA = item_surface_area

        self.FA = field_area
        self.SPAR = sparcity

        try:
            self.logSize = log2(Size)
        except:
            self.logSize = None

        try:
            self.logSpacing = log2(Spacing)
        except:
            self.logSpacing = log2(Spacing)



    def _calc(self):

        if self.TSA is None:
            self.TSA= self.N * self.ISA
        if self.TSA is None:
            self.TSA = 2**(self.logSize - log2(self.ISA))
        if self.ISA is None:
            self.ISA = float(self.TSA)/ self.N
        if self.ISA is None:
            self.ISA = 2**(self.logSize - log2(self.TSA))
        if self.FA is None:




