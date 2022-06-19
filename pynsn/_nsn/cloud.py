"""
Object Cloud
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from hashlib import md5
import json
from random import random

import numpy as np
from .._lib import misc, geometry
from .visual_features import VisualFeatures

class _Cloud(object):
    """Numpy Position lists with attributes for optimized for numpy calculations

    Abstract class for implementation of dot and rect
    """

    def __init__(self, xy=None, attributes=None):

        self._xy = np.array([])
        self._attributes = np.array([])
        self._features = VisualFeatures(self)

        if xy is not None:
            self._append_xy_attribute(xy=xy, attributes=attributes)

    def __str__(self):
        return self._features.get_features_text(extended_format=True)

    @property
    def xy(self):
        return self._xy

    @property
    def xy_rounded_integer(self):
        """rounded to integer"""
        return np.array(np.round(self._xy))

    @property
    def attributes (self):
        return self._attributes

    @property
    def features(self):
        return self._features

    def _append_xy_attribute(self, xy, attributes=None):
        """returns number of added rows"""
        xy = misc.numpy_array_2d(xy)
        if not isinstance(attributes, (tuple, list)):
            attributes = [attributes] * xy.shape[0]

        if len(attributes) != xy.shape[0]:
            raise RuntimeError(u"Bad shaped data: " + u"attributes have not "
                                                      u"the same length as the coordinates")

        self._attributes = np.append(self._attributes, attributes)
        if len(self._xy) == 0:
            empty = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
            self._xy = np.append(empty, xy, axis=0)
        else:
            self._xy = np.append(self._xy, xy, axis=0)
        self._features.reset()
        return xy.shape[0]

    def set_attributes(self, attributes):
        """Set all attributes

        Parameter
        ---------
        attributes:  attribute (string) or list of attributes

        """

        if isinstance(attributes, (list, tuple)):
            if len(attributes) != self._features.numerosity:
                raise ValueError("Length of attribute list does not match the " +\
                                 "size of the dot array.")
            self._attributes = np.array(attributes)
        else:
            self._attributes = np.array([attributes] * self._features.numerosity)

    @property
    def center_of_outer_positions(self): #FIXME it is the center of positions
        minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def surface_areas(self):
        # all zero
        return np.zeros(self._xy.shape[0])

    @property
    def perimeter(self):
        # all zero
        return np.zeros(self._xy.shape[0])

    @property
    def hash(self):
        """md5_hash of positions and perimeter"""

        m = md5()
        m.update(self._xy.tobytes())  # to byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        m.update(self.perimeter.tobytes())
        return m.hexdigest()

    def as_dict(self):
        """
        """
        d = {"xy": self._xy.tolist()}
        if misc.is_all_equal(self._attributes):
            d.update({"attributes": self._attributes[0]})
        else:
            d.update({"attributes": self._attributes.tolist()})
        return d

    def read_from_dict(self, dict):
        """read dot array from dict"""
        self._xy = np.array(dict["xy"])
        if not isinstance(dict["attributes"], (list, tuple)):
            att = [dict["attributes"]] * self._features.numerosity
        else:
            att = dict["attributes"]
        self._attributes = np.array(att)
        self._features.reset()

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
        self._features.reset()

    def clear(self):
        self._xy = np.array([[]])
        self._attributes = np.array([])
        self._features.reset()

    def delete(self, index):
        self._xy = np.delete(self._xy, index, axis=0)
        self._attributes = np.delete(self._attributes, index)
        self._features.reset()

    def __jitter_identical_positions(self, jitter_size=0.1):
        """jitters points with identical position"""

        for idx, ref_object in enumerate(self._xy):
            identical = np.where(np.all(np.equal(self._xy, ref_object), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != idx:
                        self._xy[x, :] = self._xy[x, :] - geometry.polar2cartesian(
                            [[jitter_size, random() * 2 * np.pi]])[0]

    def _radial_replacement_from_reference_dots(self, ref_pos_id,
                                                neighbour_ids, replacement_size):
        """remove neighbouring position radially from reference position
        helper function, typically used for realign
        """

        # check if there is an identical position and jitter to avoid fully overlapping positions
        if np.sum(np.all(self._xy[neighbour_ids,] == self._xy[ref_pos_id, :],
                       axis=1)) > 0:
            self.__jitter_identical_positions()

        # relative polar positions to reference_dot
        tmp_polar = geometry.cartesian2polar(self._xy[neighbour_ids, :] - self._xy[ref_pos_id, :])
        tmp_polar[:, 0] = 0.000000001 + replacement_size # determine movement size
        xy = geometry.polar2cartesian(tmp_polar)
        self._xy[neighbour_ids, :] = np.array([self._xy[neighbour_ids, 0] + xy[:, 0],
                                               self._xy[neighbour_ids, 1] + xy[:, 1]]).T


class _RestrictedCloud(_Cloud):

    def __init__(self, target_array_radius, minimum_gap,
                 xy=None, attributes=None):
        """Just a Cloud with additional variables target_array_radius, minimum_gap
         used as parent for Object Arrays"""
        super().__init__(xy=xy, attributes=attributes)

        self.target_array_radius = target_array_radius
        self.minimum_gap = minimum_gap
