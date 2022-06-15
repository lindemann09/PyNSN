"""
Dot Array
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

# TODO: How to deal with rounding? Is saving to precises? Suggestion:
#  introduction precision parameter that is used by as_dict and get_csv and
#  hash

class _DotCloud(object):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    def __init__(self, xy=None, diameters=None):

        self._xy = np.array([])
        self._diameters = np.array([])
        self._match = FeatureMatcher(self)
        self._features = VisualFeatures(self)

        if (xy, diameters) != (None, None):
            self.append(xy, diameters)

    def __str__(self):
        return self.features.get_features_text(extended_format=True)

    @property
    def xy(self):
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
    def center_of_outer_positions(self):
        minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def surface_areas(self):
        # a = pi r**2 = pi d**2 / 4
        return np.pi * (self._diameters ** 2) / 4.0

    @property
    def perimeter(self):
        return np.pi * self._diameters

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
        self._diameters = np.round(self._diameters, decimals=decimals)
        if decimals==0:
            self._xy = self._xy.astype(int_type)
            self._diameters = self._diameters.astype(int_type)

    def json(self, indent=None, include_hash=False):
        """"""
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

        with open(json_file_name, 'r') as fl:
            dict = json.load(fl)
        self.read_from_dict(dict)

    def as_dict(self):
        """
        """
        return {"xy": self._xy.tolist(),
                "diameters": self._diameters.tolist()}

    def read_from_dict(self, dict):
        """read Dot collection from dict"""
        self._xy = np.array(dict["xy"])
        self._diameters = np.array(dict["diameters"])
        self.features.reset()

    @property
    def diameters(self):
        return self._diameters

    def append(self, xy, item_diameters):
        """append dots using numpy array"""

        # ensure numpy array
        item_diameters = misc.numpy_vector(item_diameters)
        # ensure xy is a 2d array
        xy = np.array(xy)
        if xy.ndim == 1 and len(xy) == 2:
            xy = xy.reshape((1, 2))
        if xy.ndim != 2:
            raise RuntimeError("Bad shaped data: xy must be pair of xy-values or a list of xy-values")

        if xy.shape[0] != len(item_diameters):
            raise RuntimeError("Bad shaped data: " + u"xy has not the same length as item_diameters")

        if len(self._xy) == 0:
            self._xy = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
        self._xy = np.append(self._xy, xy, axis=0)
        self._diameters = np.append(self._diameters, item_diameters)
        self.features.reset()

    def clear(self):
        self._xy = np.array([[]])
        self._diameters = np.array([])
        self.features.reset()

    def delete(self, index):
        self._xy = np.delete(self._xy, index, axis=0)
        self._diameters = np.delete(self._diameters, index)
        self.features.reset()

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """
        if indices is None:
            indices = list(range(self.features.numerosity))

        return _DotCloud(xy=self._xy[indices, :].copy(),
                         diameters=self._diameters[indices].copy())

    def distances(self, xy, diameter):
        """Distances toward a single point (xy, diameter) """
        if len(self._xy) == 0:
            return np.array([])
        else:
            return np.hypot(self._xy[:, 0] - xy[0], self._xy[:, 1] - xy[1]) - \
                   ((self._diameters + diameter) / 2.0)

    @property
    def center_of_mass(self):
        weighted_sum = np.sum(self._xy * self._diameters[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self._diameters)


    def center_array(self):
        self._xy = self._xy - self.center_of_outer_positions
        self.features.reset()

    def append_dot(self, dot):
        assert isinstance(dot, Dot)
        self.append(xy=[dot.x, dot.y],
                    item_diameters=dot.diameter)

    def get_dots(self, indices=None):
        """returns all dots

        indices int or list of ints
         filtering possible:
         if diameter is defined it returns only dots a particular
         diameter

         """

        if indices is not None:
            indices = range(self.features.numerosity)
        try:
            indices = list(indices)  # check if iterable
        except:
            indices = [indices]

        rtn = []
        for xy, dia in zip(self._xy[indices, :], self._diameters[indices]):
            rtn.append(Dot(x=xy[0], y=xy[1], diameter=dia))
        return rtn



class _GenericDotArray(_DotCloud):
    """Generic Dot array is restricted to a certain area, but it has not attributes"""


    def __init__(self, target_array_radius, minimum_gap):
        """Dot array is restricted to a certain area, it has a target area
        and a minimum gap.
        This features allow find shuffling free position and matching
        features.

        target_array_radius or dot_array_file needs to be define."""

        super().__init__()
        self.target_array_radius = target_array_radius
        self.minimum_gap = minimum_gap

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        rtn = super().copy(indices=indices)
        rtn.target_array_radius = self.target_array_radius
        rtn.minimum_gap = self.minimum_gap
        return rtn

    def _jitter_identical_positions(self, jitter_size=0.1):
        """jitters points with identical position"""

        for idx, ref_dot in enumerate(self._xy):
            identical = np.where(np.all(np.equal(self._xy, ref_dot), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != idx:
                        self._xy[x, :] = self._xy[x, :] - geometry.polar2cartesian(
                            [[jitter_size, random.random() * 2 * np.pi]])[0]

    def _remove_overlap_for_dot(self, dot_id, minimum_gap):
        """remove overlap for one point
        helper function, please use realign
        """

        dist = self.distances(self._xy[dot_id, :], self._diameters[dot_id])

        shift_required = False
        idx = np.where(dist < minimum_gap)[0].tolist()  # overlapping dot ids
        if len(idx) > 1:
            idx.remove(dot_id)  # don't move yourself
            if np.sum(
                    np.all(self._xy[idx,] == self._xy[dot_id, :],
                           axis=1)) > 0:  # check if there is an identical position
                self._jitter_identical_positions()

            tmp_polar = geometry.cartesian2polar(self._xy[idx, :] - self._xy[dot_id, :])
            tmp_polar[:, 0] = 0.000000001 + minimum_gap - dist[idx]  # determine movement size
            xy = geometry.polar2cartesian(tmp_polar)
            self._xy[idx, :] = np.array([self._xy[idx, 0] + xy[:, 0], self._xy[idx, 1] + xy[:, 1]]).T
            shift_required = True

        return shift_required

    def remove_overlap_from_inner_to_outer(self, minimum_gap):

        shift_required = False
        # from inner to outer remove overlaps
        for i in np.argsort(geometry.cartesian2polar(self._xy, radii_only=True)):
            if self._remove_overlap_for_dot(dot_id=i, minimum_gap=minimum_gap):
                shift_required = True

        if shift_required:
            self.features.reset()

        return shift_required

    def realign(self):
        """Realigns the dots in order to remove all dots overlaps and dots
        outside the target area.

        If two dots overlap, the dots that is further apart from the array
        center will be moved opposite to the direction of the other dot until
        there is no overlap (note: minimun_gap parameter). If two dots have
        exactly the same position the same position one is minimally shifted
        in a random direction.

        Note: Realignming might change the field area! Match Space parameter after
        realignment.

        """

        error = False

        shift_required = self.remove_overlap_from_inner_to_outer(minimum_gap=self.minimum_gap)

        # sqeeze in points that pop out of the image area radius
        cnt = 0
        while True:
            radii = geometry.cartesian2polar(self._xy, radii_only=True)
            too_far = np.where((radii + self._diameters // 2) > self.target_array_radius)[0]  # find outlier
            if len(too_far) > 0:

                # squeeze in outlier
                polar = geometry.cartesian2polar([self._xy[too_far[0], :]])[0]
                polar[0] = self.target_array_radius - self._diameters[
                    too_far[0]] // 2 - 0.000000001  # new radius #todo check if 0.00001 required
                new_xy = geometry.polar2cartesian([polar])[0]
                self._xy[too_far[0], :] = new_xy

                # remove overlaps centered around new outlier position
                self._xy = self._xy - new_xy
                # remove all overlaps (inner to outer, i.e. starting with outlier)
                self.remove_overlap_from_inner_to_outer(minimum_gap=self.minimum_gap)
                # new pos for outlier
                self._xy = self._xy + new_xy  # move back to old position
                shift_required = True
            else:
                break  # end while loop

            cnt += 1
            if cnt > 20:
                error = True
                break

        if error:
            return False, u"Can't find solution when removing outlier (n=" + \
                   str(self.features.numerosity) + ")"

        self.features.reset()
        if not shift_required:
            return True, ""
        else:
            return self.realign()  # recursion

    def random_free_dot_position(self, dot_diameter,
                                 allow_overlapping=False,
                                 prefer_inside_field_area=False,
                                 squared_array = False,
                                 min_distance_area_boarder = 0,
                                 occupied_space=None):
        """returns a available random xy position

        raise exception if not found
        occupied space: see generator generate
        """

        try_out_inside_convex_hull = 1000

        if prefer_inside_field_area:
            delaunay = spatial.Delaunay(self.features.convex_hull.xy)
        else:
            delaunay = None
        cnt = 0

        target_radius = (self.target_array_radius -
                         min_distance_area_boarder) - (dot_diameter / 2.0)
        while True:
            cnt += 1
            ##  polar method seems to produce centeral clustering
            #  proposal_polar =  np.array([random.random(), random.random()]) *
            #                      (target_radius, TWO_PI)
            #proposal_xy = misc.polar2cartesian([proposal_polar])[0]
            # note: np.random produceing identical numbers under multiprocessing

            proposal_xy = np.array([random.random(), random.random()]) \
                          * 2 * self.target_array_radius - self.target_array_radius

            bad_position = False
            if not squared_array:
                bad_position =  target_radius <= \
                                np.hypot(proposal_xy[0], proposal_xy[1])

            if not bad_position and prefer_inside_field_area and cnt < \
                    try_out_inside_convex_hull:
                bad_position = delaunay.find_simplex(proposal_xy) < 0

            if not bad_position and not allow_overlapping:
                # find bad_positions
                dist = self.distances(proposal_xy, dot_diameter)
                if occupied_space:
                    dist = np.append(dist, occupied_space.distances(proposal_xy, dot_diameter))
                idx = np.where(dist < self.minimum_gap)[0]  # overlapping dot ids
                bad_position = len(idx) > 0

            if not bad_position:
                return proposal_xy
            elif cnt > 3000:
                raise RuntimeError(u"Can't find a free position") #FIXME

    def shuffle_all_positions(self, allow_overlapping=False,
                              min_distance_area_boarder=0):
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        new_xy = None
        for d in self.diameters:
            try:
                xy = self.random_free_dot_position(d,
                                                   allow_overlapping=allow_overlapping,
                                                   min_distance_area_boarder=min_distance_area_boarder)
            except:
                raise RuntimeError("Can't shuffle dot array. No free positions.")

            if new_xy is None:
                new_xy = np.array([xy])
            else:
                new_xy = np.append(new_xy, [xy], axis=0)

        self._xy = new_xy
        self.features.reset()

    def number_deviant(self, change_numerosity, prefer_keeping_field_area=False):
        """number deviant
        """

        TRY_OUT = 100
        # make a copy for the deviant
        deviant = self.copy()
        if self.features.numerosity + change_numerosity <= 0:
            deviant.clear()
        else:
            # add or remove random dots
            for _ in range(abs(change_numerosity)):
                if prefer_keeping_field_area:
                    ch = deviant.features.convex_hull.indices
                else:
                    ch = []
                for x in range(TRY_OUT):##
                    rnd = random.randint(0, deviant.features.numerosity - 1) # do not use np.random
                    if rnd not in ch or change_numerosity > 0:
                        break

                if change_numerosity < 0:
                    # remove dots
                    deviant.delete(rnd)
                else:
                    # copy a random dot
                    rnd_dot = self.get_dots([rnd])[0]
                    try:
                        rnd_dot.xy = deviant.random_free_dot_position(
                            dot_diameter=rnd_dot.diameter,
                            prefer_inside_field_area=prefer_keeping_field_area)
                    except:
                        # no free position
                        raise RuntimeError("Can't make the deviant. No free position")
                    deviant.append_dot(rnd_dot)

        return deviant

    def as_dict(self):
        d = super().as_dict()
        d.update({"minimum_gap": self.minimum_gap,
             "target_array_radius": self.target_array_radius})
        return d

    def read_from_dict(self, dict):
        super().read_from_dict(dict)
        self.minimum_gap = dict["minimum_gap"]
        self.target_array_radius = dict["target_array_radius"]


class DotArray(_GenericDotArray):
    """Dot array with attributes
    """

    def __init__(self, target_array_radius, minimum_gap):
        """ Dot array adds attributes, such as colour to a SimpleDotArray."""

        super().__init__(target_array_radius=target_array_radius,
                         minimum_gap=minimum_gap)
        self._attributes = np.array([])

    @property
    def attributes(self):
        return self._attributes

    def set_attributes(self, attributes):
        """Set all attributes

        Parameter
        ---------
        attributes:  attribute (string) or list of attributes

        """

        if isinstance(attributes, (list, tuple)):
            if len(attributes) != self.features.numerosity:
                raise ValueError("Length of attribute list does not match the " +\
                                 "size of the dot array.")
            self._attributes = np.array(attributes)
        else:
            self._attributes = np.array([attributes] * self.features.numerosity)

    def clear(self):
        super().clear()
        self._attributes = np.array([])

    def append(self, xy, item_diameters, attributes=None):
        """append dots with attribvute"""

        item_diameters = misc.numpy_vector(item_diameters)
        super().append(xy=xy, item_diameters=item_diameters)

        if not isinstance(attributes, (tuple, list)):
            attributes = [attributes] * len(item_diameters)

        if len(attributes) != len(item_diameters):
            raise RuntimeError(u"Bad shaped data: " + u"attributes have not "
                                                      u"the same length as diameter")
        self._attributes = np.append(self._attributes, attributes)


    def delete(self, index):
        super().delete(index)
        self._attributes = np.delete(self._attributes, index)

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """
        if indices is None:
            indices = list(range(self.features.numerosity))

        rtn = DotArray(target_array_radius=self.target_array_radius,
                        minimum_gap=self.minimum_gap)

        for d in self.get_dots(indices=indices):
            rtn.append_dot(d)
        return rtn

    def append_dot(self, dot):
        assert isinstance(dot, Dot)
        self.append(xy=[dot.x, dot.y],
                    item_diameters=dot.diameter,
                    attributes=dot.attribute)

    def get_dots(self, indices=None, diameter=None, item_attributes=None):
        """returns all dots
        indices : int or list of ins

         filtering possible:
         if diameter/colour/picture is defined it returns only dots a
         particular diameter/colour/picture

         """
        rtn = []
        i = -1
        if indices is not None:
            try:
                indices = list(indices)  # check if iterable
            except:
                indices = [indices]

        for xy, dia, att in zip(self._xy, self._diameters, self._attributes):
            i += 1
            if (indices is not None and i not in indices) or \
                    (diameter is not None and dia != diameter) or \
                    (item_attributes is not None and att !=  item_attributes):
                continue

            rtn.append(Dot(x=xy[0], y=xy[1], diameter=dia, attribute=att))
        return rtn

    def join(self, dot_array):
        """add another dot arrays"""
        self.append(xy=dot_array._xy, item_diameters=dot_array._diameters,
                    attributes=dot_array.attributes)

    def split_array_by_attributes(self):
        """returns a list of arrays
        each array contains all dots of with particular colour"""
        rtn = []
        att = self.attributes
        att[np.where(att==None )] = "None"

        for c in np.unique(att):
            if c is not None:
                da = DotArray(target_array_radius=self.target_array_radius,
                              minimum_gap=self.minimum_gap)
                for d in self.get_dots(item_attributes =c):
                    da.append_dot(d)
                rtn.append(da)
        return rtn



    def get_csv(self, variable_names=True,
                hash_column=True, num_idx_column=True,
                attribute_column=False):  # todo print features
        """Return the dot array as csv text

        Parameter
        ---------
        variable_names : bool, optional
            if True variable name will be printed in the first line

        """

        rtn = ""
        if variable_names:
            if hash_column:
                rtn += u"hash,"
            if num_idx_column:
                rtn += u"num_id,"
            rtn += u"x,y,diameter"
            if attribute_column:
                rtn += u",attribute"
            rtn += u"\n"

        obj_id = self.hash
        for cnt in range(len(self._xy)):
            if hash_column:
                rtn += "{0}, ".format(obj_id)
            if num_idx_column:
                rtn += "{},".format(self.features.numerosity)
            #rtn += num_format % self._xy[cnt, 0] + "," + num_format %
                # self._xy[cnt, 1] + "," + \
            #       num_format % self._diameters[cnt]
            rtn += "{},{},{}".format(self._xy[cnt, 0], self._xy[cnt, 1],
                                     self._diameters[cnt])
            if attribute_column:
                rtn += ", {}".format(self.attributes[cnt])
            rtn += "\n"
        return rtn

    def as_dict(self):
        d = super().as_dict()
        if misc.is_all_equal(self._attributes):
            d.update({"attributes": self._attributes[0]})
        else:
            d.update({"attributes": self._attributes.tolist()})
        return d

    def read_from_dict(self, d):
        super().read_from_dict(d)
        if not isinstance(d["attributes"], (list, tuple)):
            att = [d["attributes"]] * self.features.numerosity
        else:
            att = d["attributes"]

        self._attributes = np.array(att)
