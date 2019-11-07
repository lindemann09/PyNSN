"""
Dot Array
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import numpy as np

from ._colour import Colour

class ItemAttributes(object):

    def __init__(self, colour=None, picture=None):
        self.colour = Colour(colour)
        self.picture = picture

    def as_dict(self):
        return {"color" : str(self.colour),
                "picture" : self.picture}


class ItemAttributesList(object):

    def __init__(self):
        """If all ItemAttributes should None (but not empty) init with ItemAttributes(colours=[None]),
        otherwise the ItemAttributes will be empty"""

        self.clear()

    def clear(self):
        self.colours = np.array([])
        self.pictures = np.array([])

    def __iter__(self):
        return map(lambda x: ItemAttributes(colour=x[0], picture=x[1]),
                   zip(self.colours, self.pictures))

    def __getitem__(self, key):
        return ItemAttributes(colour=self.colours[key],
                              picture=self.pictures[key])

    def as_dict(self):

        rtn = {"colours": list(map(lambda x: x._colour, self.colours))}

        if not np.all(self.pictures == None):
            rtn.update({"pictures": self.pictures.tolist()})
        return rtn

    def read_from_dict(self, dict):
        col = list(map(lambda x: Colour(x), dict["colours"]))
        self.colours = np.array(col)
        try:
            self.pictures = np.array(dict["pictures"])
        except:
            self.pictures = np.array([None] * len(self.colours))

    @property
    def length(self):
        return len(self.colours)

    def delete(self, index):
        self.colours = np.delete(self.colours, index)
        self.pictures = np.delete(self.pictures, index)

    def repeat(self, repeats):
        """:return repeated attributes, equivalent to numpy.repeat"""

        rtn = ItemAttributesList()
        rtn.colours = self.colours.repeat(repeats)
        rtn.pictures = self.pictures.repeat(repeats)
        return rtn

    def copy(self):
        rtn = ItemAttributesList()
        rtn.colours = self.colours.copy()
        rtn.pictures = self.pictures.copy()
        return rtn

    def append(self, attributes):
        """attributes: ItemAttributes, ItemAttributesList or list of ItemAttributes
        """

        if isinstance(attributes, ItemAttributes):
            self._append_lists(colours=[attributes.colour],
                               pictures=[attributes.picture])

        elif isinstance(attributes, (list, tuple)):
            # list of ItemAttributes
            for att in attributes:
                self.append(att)

        elif isinstance(attributes, ItemAttributesList):
            self._append_lists(colours=attributes.colours,
                               pictures=attributes.pictures)
        else:
            raise ValueError("arrtibutes need to be ItemAttributes, list "
                             "of ItemAttributes or ItemAttributesList")



    def _append_lists(self, colours, pictures):
        # lists of single attributes with equal length

        if (len(colours) != len(pictures)):
            raise RuntimeError(u"Bad shaped data: attributes have not the same length.")

        self.pictures = np.append(self.pictures, pictures)  # TODO check picture exist
        self.colours = np.append(self.colours, list(map(lambda x:Colour(x),
                                                    colours)))

    def change(self, colour=None, picture=None, indices=None):
        """allows using color names and rgb array, since it
        converts colour """

        if isinstance(indices, int):
            indices = [indices]
        elif indices is None:
            indices = range(len(self.colours))

        if colour is not None:
            self.colours[indices] = Colour(colour)
        if picture is not None:
            self.pictures[indices] = picture


