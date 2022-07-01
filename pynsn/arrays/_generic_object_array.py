"""
Object Cloud
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from hashlib import md5
import json

import numpy as np
from scipy import spatial
from .._lib import misc
from ..visual_properties._props import ArrayProperties
from ..visual_properties import fit

class GenericObjectArray(object):

    def __init__(self, target_area_radius,
                 min_dist_between = 2,
                 min_dist_area_boarder = 1,
                 xy=None,
                 attributes=None):
        """Numpy Position lists with attributes for optimized for numpy calculations

        Abstract class for implementation of dot and rect
        """
        self.target_area_radius = target_area_radius
        self.min_dist_between = min_dist_between
        self.min_dist_area_boarder = min_dist_area_boarder

        self._xy = np.array([])
        self._attributes = np.array([])
        self._properties = ArrayProperties(self)

        if xy is not None:
            self._append_xy_attribute(xy=xy, attributes=attributes)


    def __str__(self):
        return self._properties.as_text(extended_format=True)

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
    def properties(self):
        return self._properties

    def _append_xy_attribute(self, xy, attributes=None):
        """returns number of added rows"""
        xy = misc.numpy_array_2d(xy)
        if not isinstance(attributes, (tuple, list)):
            attributes = [attributes] * xy.shape[0]

        if len(attributes) != xy.shape[0]:
            raise ValueError(u"Bad shaped data: " + u"attributes have not "
                                                      u"the same length as the coordinates")

        self._attributes = np.append(self._attributes, attributes)
        if len(self._xy) == 0:
            empty = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
            self._xy = np.append(empty, xy, axis=0)
        else:
            self._xy = np.append(self._xy, xy, axis=0)
        self._properties.reset()
        return xy.shape[0]

    def set_attributes(self, attributes):
        """Set all attributes

        Parameter
        ---------
        attributes:  attribute (string) or list of attributes

        """

        if isinstance(attributes, (list, tuple)):
            if len(attributes) != self._properties.numerosity:
                raise ValueError("Length of attribute list does not adapt the " +\
                                 "size of the dot array.")
            self._attributes = np.array(attributes)
        else:
            self._attributes = np.array([attributes] * self._properties.numerosity)

    def center_of_positions(self):
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

    def get_number_deviant(self, change_numerosity,
                           keeping_field_area=False):
        """number deviant
        """
        object_array = self.copy()
        new_num = self.properties.numerosity + change_numerosity
        fit.numerosity(object_array, value=new_num, keeping_field_area=False)
        return object_array

    def as_dict(self):
        """
        """
        d = {"xy": self._xy.tolist()}
        if len(self._attributes) >0 and misc.is_all_equal(self._attributes):
            d.update({"attributes": self._attributes[0]})
        else:
            d.update({"attributes": self._attributes.tolist()})
        return d

    def read_from_dict(self, dict):
        """read dot array from dict"""
        self._xy = np.array(dict["xy"])
        if not isinstance(dict["attributes"], (list, tuple)):
            att = [dict["attributes"]] * self._properties.numerosity
        else:
            att = dict["attributes"]
        self._attributes = np.array(att)
        self._properties.reset()

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
        self._xy = self._xy - self.center_of_positions()
        self._properties.reset()

    def clear(self):
        self._xy = np.array([[]])
        self._attributes = np.array([])
        self._properties.reset()

    def delete(self, index):
        self._xy = np.delete(self._xy, index, axis=0)
        self._attributes = np.delete(self._attributes, index)
        self._properties.reset()

    def specifications_dict(self):
        return {"target_area_radius": self.target_area_radius,
                "min_dist_between": self.min_dist_between,
                "min_dist_area_boarder": self.min_dist_area_boarder}

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if indices is None:
            indices = list(range(self._properties.numerosity))

        return GenericObjectArray(target_area_radius=self.target_area_radius,
                        min_dist_between=self.min_dist_between,
                        min_dist_area_boarder = self.min_dist_area_boarder,
                        xy=self._xy[indices, :].copy(),
                        attributes=self._attributes[indices].copy())

    def get(self):
        raise NotImplementedError()

    def center_of_mass(self):
        weighted_sum = np.sum(self._xy * self.perimeter[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self.perimeter)

    def distances(self, ref_object):
        # override this method
        raise NotImplementedError()

    def distances_matrix(self, between_positions=False, overlap_is_zero=False):
        """between position ignores the dot size"""
        if between_positions:
            return spatial.distance.cdist(self._xy, self._xy)
        # matrix with all distance between all points
        dist = np.array([self.distances(d) for d in self.get()])
        if overlap_is_zero:
            dist[dist<0] = 0
        return dist

