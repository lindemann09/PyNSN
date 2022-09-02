from __future__ import annotations
from zlib import DEF_BUF_SIZE

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from os import path
from typing import Any, Iterator, Union

from numpy.typing import ArrayLike, NDArray

from .rectangle import Rectangle


class Picture(Rectangle):
    ATTR_PREFIX = "p:"

    def __init__(self,
                 xy: ArrayLike = (0, 0),
                 size: ArrayLike = (0, 0),
                 filename: str = ""):
        """Initialize a Picture

        Rectangle can also consist of a picture

        Parameters
        ----------
        xy : tuple
            tuple of two numeric (default=(0, 0))
        size : tuple
            tuple of two numeric (default=(0, 0))
        attribute : attribute (string) or PictureFile
        """


        super().__init__(xy=xy, size=size,
                         attribute=Picture.ATTR_PREFIX + filename)

    @property
    def filename(self):
        return self.attribute[len(Picture.ATTR_PREFIX):]

    def check_file_exists(self):
        """Checks if the file exists
        """
        return path.isfile(self.filename)

    @staticmethod
    def is_picture_attribute(txt: str) -> bool:
        """Check if text string is a picture attribute

        Args:
            txt: string to be checked
        """
        return isinstance(txt, str) and \
                txt.startswith(Picture.ATTR_PREFIX)

