"""
Object Cloud
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import random
from hashlib import md5
import json

import numpy as np
from scipy import spatial

from ..lib import misc, geometry
from pynsn.dot_array.shape import Dot
from .visual_features import VisualFeatures
from .match import FeatureMatcher


def np_array_2d(two_d_data): #FIXME move to misc
    """ensures well shaped to numpy array"""
    rtn = np.array(two_d_data)
    if rtn.ndim == 1 and len(rtn) == 2:
        rtn = rtn.reshape((1, 2))
    if rtn.ndim != 2:
        raise RuntimeError("Bad shaped data: xy must be pair of xy-values or a list of xy-values")
    return rtn


class _ObjectCloud(object):
    """Numpy Position lists for optimized for numpy calculations

    Abstract class for implementation of dot and rect
    """

    def __init__(self):
        self._xy = np.array([])
        self._match = FeatureMatcher(self)
        self._features = VisualFeatures(self)

    def __str__(self):
        return self.features.get_features_text(extended_format=True)

    @property
    def xy(self):
        if len(self._xy) == 0:
            return np.array([]).reshape((0, 2))  # ensure good shape of self.xy
        else:
            return self._xy

    @property
    def match(self):
        return self._match

    @property
    def features(self):
        return self._features

    @property
    def xy_rounded_integer(self):
        """rounded to integer"""
        return np.array(np.round(self._xy))

    @property
    def center_of_outer_positions(self): #FIXME it is the center of positions
        minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def surface_areas(self):
        return None #FIXME abstract

    @property
    def hash(self):
        """md5_hash of position, diameter"""

        m = md5()
        m.update(self._xy.tobytes())  # to byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        m.update(self.surface_areas.tobytes())
        return m.hexdigest()

    def round(self, decimals=0, int_type=np.int64):
        """Round values of the array."""

        if decimals is None:
            return
        self._xy = np.round(self._xy, decimals=decimals)
        if decimals == 0:
            self._xy = self._xy.astype(int_type)

    def as_dict(self):
        """
        """
        return {"xy": self._xy.tolist()}

    def read_from_dict(self, dict):
        """read Dot collection from dict"""
        self._xy = np.array(dict["xy"])
        self.features.reset()

    def json(self, indent=None, include_hash=False):
        """"""
        # override and extend as_dict not this function

        d = self.as_dict()
        if include_hash:
            d.update({"hash": self.hash})
        if not indent:
            indent = None
        return json.dumps(d, indent=indent)

    def save(self, json_file_name, indent=None, include_hash=False):
        """"""
        with open(json_file_name, 'w') as fl:
            fl.write(self.json(indent=indent, include_hash=include_hash))

    def load(self, json_file_name):
        # override and extend read_from_dict not this function
        with open(json_file_name, 'r') as fl:
            dict = json.load(fl)
        self.read_from_dict(dict)

    def center_array(self):
        self._xy = self._xy - self.center_of_outer_positions
        self.features.reset()

    def clear(self):
        self._xy = np.array([[]])
        self.features.reset()

    def delete(self, index):
        self._xy = np.delete(self._xy, index, axis=0)
        self.features.reset()

