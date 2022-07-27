"""

"""
from __future__ import annotations

from os import path


class PictureFile(object):
    ATTR_PREFIX = "p:"

    def __init__(self, filename: str) -> None:
        """Picture File

        Args:
            filename: filename or attribute string of a filename("p:<filename>")
        """
        if PictureFile.is_picture_attribute(filename):
            self._filename = filename[len(PictureFile.ATTR_PREFIX):]
        else:
            self._filename = filename

    @property
    def filename(self) -> str:
        """Filename (str)
        """
        return self._filename

    @property
    def attribute(self):
        """Attribute string of the filename

        which is 'p:<filename>', see ``PictureFile.ATTR_PREFIX``
        """
        return "f{PictureFile.PICTURE_PREFIX}f{self._filename}"

    @staticmethod
    def is_picture_attribute(txt: str) -> bool:
        """Check if text string is a picture attribute

        Args:
            txt: string to be checked
        """
        return isinstance(txt, str) and \
                txt.startswith(PictureFile.ATTR_PREFIX)

    def check_file_exists(self):
        """Checks if the file exists
        """
        return path.isfile(self.filename)
