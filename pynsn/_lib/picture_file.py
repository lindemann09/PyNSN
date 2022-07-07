from os import path

class PictureFile(object):
    PICTURE_PREFIX = "p:"

    def __init__(self, filename):
        self._attribute = "{}{}".format(PictureFile.PICTURE_PREFIX, filename)

    @property
    def filename(self):
        return PictureFile.check_attribute(self._attribute)

    @property
    def attribute(self):
        return self._attribute

    @classmethod
    def check_attribute(cls, txt):
        """returns filename if `txt` (str) represent a picture attribute
        else returns None
        """
        if isinstance(txt, str) and \
                txt.startswith(PictureFile.PICTURE_PREFIX):
            return txt[len(PictureFile.PICTURE_PREFIX):]
        else:
            return None

    def check_file_exists(self):
        return path.isfile(self.filename)
