__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod

from ..coordinate import Coordinate
from ...image._colour import Colour
from .picture_file import PictureFile


# todo typing


class ABCShape(Coordinate, metaclass=ABCMeta):

    def __init__(self, xy, attribute):
        Coordinate.__init__(self, x=xy[0], y=xy[1])
        self._attribute = None
        self.attribute = attribute  # call setter

    @property
    def attribute(self):
        return self._attribute

    @attribute.setter
    def attribute(self, attr):
        """set attribute 
        
        Parameters
        ----------
        attr : anything
            can be, in principle, anything.
            If Colour or PictureFile, it will convert it to their string
            representation
        """
        if isinstance(attr, Colour):
            self._attribute = attr.colour
        elif isinstance(attr, PictureFile):
            self._attribute = attr.attribute
        else:
            self._attribute = attr

    def get_attribute_object(self):
        """Class instance of the attribute, if possible

        Returns
        -------
        rtn : attribute
            If attribute represents Colour or PictureFile it returns the instance
            of the respective class otherwise the previously defined  attribute
        """

        if isinstance(self._attribute, str):
            # check is color or picture
            col = Colour(self._attribute)
            if col.colour is not None:
                return col
            else:
                fl_name = PictureFile.check_attribute(self._attribute)
                if fl_name is not None:
                    return PictureFile(fl_name)

        return self._attribute

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def distance(self, other):
        pass

    @property
    @abstractmethod
    def area(self):
        pass

    @property
    @abstractmethod
    def perimeter(self):
        pass
