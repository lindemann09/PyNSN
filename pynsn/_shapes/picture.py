__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from os import path
from typing import Union

from .._lib.geometry import Coord2DLike
from .rectangle import Rectangle


class Picture(Rectangle):
    ATTR_PREFIX = "p:"

    def __init__(
        self, filename: str = "", xy: Coord2DLike = (0, 0), size: Coord2DLike = (0, 0)
    ) -> None:
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

        super().__init__(xy=xy, size=size, attribute=Picture.ATTR_PREFIX + filename)

    def __repr__(self):
        return (
            f"Picture(xy={self.xy}, size={self.size}, " + f"filename='{self.filename}')"
        )

    @property
    def filename(self) -> str:
        return self.attribute[len(Picture.ATTR_PREFIX) :]

    def file_exists(self) -> bool:
        """Checks if the file exists"""
        return path.isfile(self.filename)

    @staticmethod
    def extract_filename(txt: str) -> Union[str, None]:
        """Check if text string is a picture attribute and returns the filename.
        Otherwise, returns None

        Args:
            txt: string to be checked
        """
        if isinstance(txt, str) and txt.startswith(Picture.ATTR_PREFIX):
            return txt[len(Picture.ATTR_PREFIX) :]
        else:
            return None
