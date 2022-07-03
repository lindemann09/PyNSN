__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from hashlib import md5
import json
import numpy as np
from .._lib import misc, geometry
from ..visual_properties._properties import ArrayProperties

class BaseArray(object):
    """Class for attributes on two dimensional space"""

    def __init__(self, target_area_radius,
                 min_dist_between=2,
                 min_dist_area_boarder=1,
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

    def _append_xy_attribute(self, xy, attributes=None):
        """returns number of added rows"""
        xy = misc.numpy_array_2d(xy)
        if not isinstance(attributes, (tuple, list)):
            attributes = [attributes] * xy.shape[0]

        if len(attributes) != xy.shape[0]:
            raise ValueError("Bad shaped data: attributes have not " +
                             "the same length as the coordinates")

        self._attributes = np.append(self._attributes, attributes)
        if len(self._xy) == 0:
            empty = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
            self._xy = np.append(empty, xy, axis=0)
        else:
            self._xy = np.append(self._xy, xy, axis=0)
        self._properties.reset()
        return xy.shape[0]

    def __str__(self):
        prop_text = self._properties.as_text(extended_format=True)
        rtn = "- {}".format(type(self).__name__)
        rtn += "\n " + prop_text[1:]  # replace "-" with " "
        return rtn

    @property
    def xy(self):
        return self._xy

    @property
    def xy_rounded_integer(self):
        """rounded to integer"""
        return np.array(np.round(self._xy))

    @property
    def attributes(self):
        return self._attributes

    @property
    def properties(self):
        return self._properties

    @property
    def surface_areas(self):
        """per definition always zero"""
        return np.array([0] * len(self._xy))

    @property
    def perimeter(self):
        """per definition always zero"""
        return np.array([0] * len(self._xy))

    def set_attributes(self, attributes):
        """Set all attributes

        Parameter
        ---------
        attributes:  attribute (string) or list of attributes

        """

        if isinstance(attributes, (list, tuple)):
            if len(attributes) != self._properties.numerosity:
                raise ValueError("Length of attribute list does not adapt the " + \
                                 "size of the dot array.")
            self._attributes = np.array(attributes)
        else:
            self._attributes = np.array([attributes] * self._properties.numerosity)

    @property
    def hash(self):
        """md5_hash of positions and perimeter"""
        m = md5()
        m.update(
            self._xy.tobytes())  # to_byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        try:
            m.update(self.perimeter.tobytes())
        except AttributeError:
            pass
        m.update(self._attributes.tobytes())
        return m.hexdigest()

    def center_of_positions(self):
        """Center of all object positions
        Notes
        -----
        if you want centre that takes into account the object size, use
        center_of_field_area
        """
        return geometry.center_of_positions(self._xy)

    def clear(self):
        self._xy = np.array([])
        self._attributes = np.array([])
        self._properties.reset()

    def delete(self, index):
        self._xy = np.delete(self._xy, index, axis=0)
        self._attributes = np.delete(self._attributes, index)
        self._properties.reset()

    def copy(self, indices=None, deepcopy=True):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._xy) == 0:
            return BaseArray(target_area_radius=self.target_area_radius,
                             min_dist_between=self.min_dist_between,
                             min_dist_area_boarder=self.min_dist_area_boarder)

        if indices is None:
            indices = list(range(len(self._xy)))

        if deepcopy:
            return BaseArray(target_area_radius=self.target_area_radius,
                             min_dist_between=self.min_dist_between,
                             min_dist_area_boarder=self.min_dist_area_boarder,
                             xy=self._xy[indices, :].copy(),
                             attributes=self._attributes[indices].copy())
        else:
            return BaseArray(target_area_radius=self.target_area_radius,
                             min_dist_between=self.min_dist_between,
                             min_dist_area_boarder=self.min_dist_area_boarder,
                             xy=self._xy[indices, :],
                             attributes=self._attributes[indices])

    def as_dict(self):
        """
        """
        d = {"type": type(self).__name__,
             "target_area_radius": self.target_area_radius,
             "min_dist_between": self.min_dist_between,
             "min_dist_area_boarder": self.min_dist_area_boarder}
        d.update({"xy": self._xy.tolist()})
        if len(self._attributes) > 0 and misc.is_all_equal(self._attributes):
            d.update({"attributes": self._attributes[0]})
        else:
            d.update({"attributes": self._attributes.tolist()})
        return d

    def read_from_dict(self, dict):
        """read dot array from dict"""
        self.target_area_radius = dict["target_area_radius"]
        self.min_dist_between = dict["min_dist_between"]
        self.min_dist_area_boarder = dict["min_dist_area_boarder"]
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

