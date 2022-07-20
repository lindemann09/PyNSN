__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import deepcopy
from abc import ABCMeta, abstractmethod
from hashlib import md5
import json
import numpy as np
from scipy import spatial

from . import misc
from . import geometry
from . import shapes
from . import rng
from .array_tools import BrownianMotion
from ..visual_properties._properties import ArrayProperties
from ..exceptions import NoSolutionError
from ..visual_properties import fit
from .. import constants


class ArrayParameter(object):

    def __init__(self, target_area_radius,
                 min_dist_between=None,
                 min_dist_area_boarder=None):
        """Numpy Position lists with attributes for optimized for numpy calculations

        Abstract class for implementation of dot and rect
        """
        self.target_area_radius = target_area_radius
        if min_dist_between is None:
            self.min_dist_between = constants.DEFAULT_MIN_DIST_BETWEEN
        else:
            self.min_dist_between = min_dist_between
        if min_dist_area_boarder is None:
            self.min_dist_area_boarder = constants.DEFAULT_MIN_DIST_AREA_BOARDER
        else:
            self.min_dist_area_boarder = min_dist_area_boarder

    def as_dict(self):
        return {"type": type(self).__name__,
                "target_area_radius": self.target_area_radius,
                "min_dist_between": self.min_dist_between,
                "min_dist_area_boarder": self.min_dist_area_boarder}


class AttributeArray(ArrayParameter):
    """Class for attributes on two dimensional space"""

    def __init__(self, target_area_radius,
                 min_dist_between=None,
                 min_dist_area_boarder=None,
                 xy=None,
                 attributes=None):
        """Numpy Position lists with attributes for optimized for numpy calculations

        Abstract class for implementation of dot and rect
        """
        super().__init__(target_area_radius=target_area_radius,
                         min_dist_between=min_dist_between,
                         min_dist_area_boarder=min_dist_area_boarder)

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
        return np.round(self._xy)

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

    def get_center_of_positions(self):
        """Center of all object positions
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
            return AttributeArray(target_area_radius=self.target_area_radius,
                                  min_dist_between=self.min_dist_between,
                                  min_dist_area_boarder=self.min_dist_area_boarder)

        if indices is None:
            indices = list(range(len(self._xy)))

        if deepcopy:
            return AttributeArray(target_area_radius=self.target_area_radius,
                                  min_dist_between=self.min_dist_between,
                                  min_dist_area_boarder=self.min_dist_area_boarder,
                                  xy=self._xy[indices, :].copy(),
                                  attributes=self._attributes[indices].copy())
        else:
            return AttributeArray(target_area_radius=self.target_area_radius,
                                  min_dist_between=self.min_dist_between,
                                  min_dist_area_boarder=self.min_dist_area_boarder,
                                  xy=self._xy[indices, :],
                                  attributes=self._attributes[indices])

    def as_dict(self):
        """
        """
        d = super().as_dict()
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


class ABCObjectArray(AttributeArray, metaclass=ABCMeta):

    @property
    @abstractmethod
    def surface_areas(self):
        pass

    @property
    @abstractmethod
    def perimeter(self):
        pass

    @abstractmethod
    def as_dict(self):
        return super().as_dict()

    @abstractmethod
    def read_from_dict(self, dict):
        return super().read_from_dict()

    @abstractmethod
    def copy(self, indices=None, deepcopy=True):
        pass

    @abstractmethod
    def iter_objects(self, indices=None):
        pass

    def get_objects(self, indices=None):
        return list(self.iter_objects(indices=indices))

    @abstractmethod
    def add(self, something):
        pass

    @abstractmethod
    def find_objects(self, attribute):
        pass

    @abstractmethod
    def get_distances(self, ref_object):
        # override this method
        pass

    @abstractmethod
    def csv(self):
        pass

    @abstractmethod
    def mod_round_values(self):
        pass

    def join(self, object_array):
        """add another object arrays"""
        self.add(object_array.iter_objects())

    def get_distances_matrix(self, between_positions=False):
        """between position ignores the dot size"""
        if between_positions:
            return spatial.distance.cdist(self._xy, self._xy)
        # matrix with all distance between all points
        dist = np.asarray([self.get_distances(d) for d in self.iter_objects()])
        return dist

    def get_overlaps(self):
        """return pairs of indices of overlapping of objects and an array of the
        amount of overlap
        """
        dist = np.triu(self.get_distances_matrix(between_positions=False), k=1)
        overlap = np.where(dist < 0)
        return np.asarray(overlap).T, np.abs(dist[overlap])

    def get_center_of_mass(self):
        weighted_sum = np.sum(self._xy * self.perimeter[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self.perimeter)

    def mod_center_array_mass(self):
        self._xy = self._xy - self.get_center_of_mass()
        self._properties.reset()

    def mod_center_field_area(self):
        cxy = geometry.center_of_positions(self.properties.convex_hull.xy)
        self._xy = self._xy - cxy
        self._properties.reset()

    def get_free_position(self, ref_object,
                          in_neighborhood=False,
                          allow_overlapping=False,
                          inside_convex_hull=False,
                          occupied_space=None):
        """returns the copy of object of at a random free position

        raise exception if not found
        occupied space: see generator generate
        """

        if isinstance(ref_object, shapes.Dot):
            object_size = ref_object.diameter / 2.0
        elif isinstance(ref_object, shapes.Rectangle):
            object_size = max(ref_object.size)
        else:
            raise NotImplementedError("Not implemented for {}".format(
                type(ref_object).__name__))
        if occupied_space is not None and \
                not isinstance(occupied_space, ABCObjectArray):  # FIXME check
            raise TypeError("Occupied_space has to be a Dot or Rectangle Array or None.")

        area_rad = self.target_area_radius - self.min_dist_area_boarder - object_size
        rtn_object = deepcopy(ref_object)  # tested deepcopy required

        if in_neighborhood:
            random_walk = BrownianMotion(ref_object.xy, delta=2)
        else:
            random_walk = None

        cnt = 0
        while True:
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError(u"Can't find a free position")
            cnt += 1

            if in_neighborhood:
                rtn_object.xy = random_walk.next()
            else:
                rtn_object.xy = rng.generator.random(size=2) * 2 * area_rad - area_rad

            # is outside area
            if isinstance(ref_object, shapes.Dot):
                is_outside = area_rad <= rtn_object.polar_radius
            else:
                # Rect: check if one edge is outside
                is_outside = False
                for e in rtn_object.iter_edges():
                    if e.polar_radius >= area_rad:
                        is_outside = True
                        break

            if not is_outside and inside_convex_hull:
                # use only those that do not change the convex hull
                tmp_array = self.copy(deepcopy=True)
                tmp_array.add([rtn_object])
                is_outside = tmp_array.properties.convex_hull != \
                               self.properties.convex_hull
            if is_outside:
                if in_neighborhood:
                    random_walk.step_back()
                continue

            if not allow_overlapping:
                # find overlap
                dist = self.get_distances(rtn_object)
                if isinstance(occupied_space, ABCObjectArray):
                    dist = np.append(dist, occupied_space.get_distances(rtn_object))
                if sum(dist < self.min_dist_between) > 0:  # at least one is overlapping
                    continue
            return rtn_object

    def mod_shuffle_positions(self, allow_overlapping=False):
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        all_objects = list(self.iter_objects())
        self.clear()
        for obj in all_objects:
            try:
                new = self.get_free_position(obj, in_neighborhood=False,
                                             allow_overlapping=allow_overlapping)
            except NoSolutionError as e:
                raise NoSolutionError("Can't shuffle dot array. No free positions found.")
            self.add([new])

    def get_outlier(self):
        """returns indices of object that stand out and array with the size
        of outstanding
        """

        xy = self.properties.convex_hull.xy
        sizes_outlying = np.hypot(xy[:, 0], xy[:, 1]) - \
                        (self.target_area_radius - self.min_dist_area_boarder)
        idx = sizes_outlying > 0
        return self.properties.convex_hull.object_indices[idx], sizes_outlying[idx]

    def get_number_deviant(self, change_numerosity, preserve_field_area=False):
        """number deviant
        """
        object_array = self.copy()
        new_num = self.properties.numerosity + change_numerosity
        fit.numerosity(object_array, value=new_num,
                       keep_convex_hull=preserve_field_area)
        return object_array

    def mod_remove_overlaps(self, keep_convex_hull=False, strict=False):
        """
        Returns
        Parameters
        ----------
        keep_convex_hull
        strict

        Returns
        -------
        rtn: boolean
            True, if field area has not been changed (in case strict=False)

        """

        warning_info = "Can't keep field area constant."
        old_fa = self.properties.field_area

        if not keep_convex_hull:
            # touch convex hull objects
            ids = list(range(len(self._xy)))
            for x in self.properties.convex_hull.object_indices_unique:
                self.mod_move_object(x, 0, (0, 0), push_other=True)
                ids.remove(x)
            # touch remaining ids
            for x in ids:
                # touch each object and push other
                self.mod_move_object(x, 0, (0, 0), push_other=True)

        else:
            overlaps = self.get_overlaps()[0]
            ch_idx = self.properties.convex_hull.object_indices_unique

            while len(overlaps):
                # do not replace convexhull objects and try to
                # take that one not on convex hull or raise error/warning
                if overlaps[0, 0] not in ch_idx:
                    idx = overlaps[0, 0]
                elif overlaps[0, 1] not in ch_idx:
                    idx = overlaps[0, 1]
                elif not strict:
                    # warning
                    idx = overlaps[0, 0]
                else:
                    raise NoSolutionError(warning_info)

                obj = next(self.iter_objects(idx))
                self.delete(idx)

                # search new pos: fist inside convex hull later outside
                found = None
                for inside_convex_hull in (True, False):
                    if strict and not inside_convex_hull:
                        continue
                    for in_neighborhood in (True, False):
                        if found is None:
                            try:
                                found = self.get_free_position(ref_object=obj,
                                                    in_neighborhood=in_neighborhood,
                                                    inside_convex_hull=inside_convex_hull)
                            except NoSolutionError as e:
                                found = None

                if found is None:
                    self.add([obj])
                    raise NoSolutionError("Can't find a solution for remove overlap")
                else:
                    self.add([found])
                    overlaps = self.get_overlaps()[0]

            self.properties.reset()

        # check convex hull change
        new_ch = self.properties.field_area
        if keep_convex_hull and old_fa != new_ch:
            if strict:
                raise NoSolutionError(warning_info)
            else:
                print("Warning: " + warning_info)

        return old_fa == new_ch

    def mod_move_object(self, object_id, distance, direction,
                        push_other=False):
        """

        Parameters
        ----------
        object_id
        distance
        direction: numeric (polar) or tuple, list or Point (cartesian)
            angle (numeric, polar angle coordinate) or a cartesian 2D coordinates indicating
            the direction towards the object should be moved

        Returns
        -------

        """

        try:
            ang = float(direction)
        except (TypeError, ValueError):
            try:
                ang = shapes.Point(x=direction[0], y=direction[1])
            except IndexError:
                ang = None
        if ang is None:
            raise TypeError("Direction has to be float or a 2D coordinate, that is, a "
                            "shapes.Point or  a tuple/list of two elements.")

        obj = next(self.iter_objects(indices=object_id))

        movement = shapes.Point()
        if isinstance(ang, float):
            movement.polar = (distance, ang)
        else:
            # "ang" isnot an angle it a point
            movement = ang - obj
            movement.polar_radius = distance

        self._xy[object_id, :] = self._xy[object_id, :] + movement.xy

        if push_other:
            # push overlapping object
            obj = next(self.iter_objects(indices=object_id))
            dist = self.get_distances(obj)
            for other_id in np.flatnonzero(dist<0):
                if other_id != object_id:
                    movement.xy = self._xy[other_id, :] - self._xy[object_id, :]
                    self.mod_move_object(other_id,
                                direction = movement.polar_angle,
                                distance = abs(dist[other_id]) + self.min_dist_between,
                                push_other=True)

        self.properties.reset()

    def mod_squeeze_to_area(self, push_other=True):
        """squeeze in target area to remove all standouts"""
        cnt = 0
        while True:
            cnt += 1
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError("Can't find a solution for squeezing")

            idx, size = self.get_outlier()
            if len(idx) == 0:
                return
            for object_id, size in zip(idx, size):
                self.mod_move_object(object_id, distance=size,
                                     direction = (0,0),
                                     push_other=push_other)

    def get_split_arrays(self):
        """returns a list of arrays
        each array contains all dots of with particular colour"""
        att = self._attributes
        att[np.where(att == None)] = "None"  # TODO check "is none"

        rtn = []
        for c in np.unique(att):
            if c is not None:
                da = self.copy(indices=0, deepcopy=False)  # fast. shallow copy with just one object
                da.clear()
                da.add(self.find_objects(attribute=c))
                rtn.append(da)
        return rtn

    def mod_realign(self, keep_convex_hull=False, strict=True):
        """

        Parameters
        ----------
        keep_convex_hull
        strict

        Returns
        -------
        field_area_unchanged, no_outlier

        """
        convex_hull_unchanged = self.mod_remove_overlaps(keep_convex_hull=keep_convex_hull,
                                                strict=strict)
        self.mod_center_field_area()

        has_no_outlier = len(self.get_outlier()[0]) == 0
        if not has_no_outlier:
            if keep_convex_hull:
                warning_info = "Can't keep field area constant."
                if strict:
                    raise NoSolutionError(warning_info)
                else:
                    print("Warning: " + warning_info)

            self.mod_squeeze_to_area(push_other=True)
            has_no_outlier = True
            convex_hull_unchanged = False

        return convex_hull_unchanged, has_no_outlier

# TODO  everywhere: file header doc and author information
