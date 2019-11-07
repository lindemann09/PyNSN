"""
Dot Array
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from copy import copy
import random
from collections import OrderedDict
from hashlib import md5
import json

import numpy as np
from scipy import spatial

from .geometry import Dot
from . import misc
from ._item_attributes import ItemAttributesList, ItemAttributes
from . import visual_features as vf


class DotCollection(object):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    def __init__(self, xy=None, diameters=None):

        self._xy = np.array([])
        self._diameters = np.array([])
        self._ch = None
        if (xy, diameters) != (None, None):
            self.append(xy, diameters)
        self.set_array_modified()

    def __str__(self):
        return self.get_features_text(extended_format=True)

    @property
    def xy(self):
        return self._xy

    @property
    def rounded_xy(self):
        """rounded to integer"""
        return np.array(np.round(self._xy), dtype=np.int32)

    @property
    def center_of_outer_positions(self):
        minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def feature_mean_item_diameter(self):
        return np.mean(self._diameters)

    @property
    def feature_total_surface_area(self):
        return np.sum(self.surface_areas)

    @property
    def convex_hull_positions(self):  # FIXME not defined for l<3
        return self._xy[self.convex_hull_indices, :]

    @property
    def convex_hull_indices(self):  # FIXME not defined for l<3
        """this convex_hull takes into account the dot diameter"""
        return self._ch.convex_hull.vertices

    @property
    def surface_areas(self):
        return np.pi * (self._diameters ** 2) / 4.0

    @property
    def perimeter(self):
        return np.pi * self._diameters

    @property
    def feature_mean_item_surface_area(self):
        return np.mean(self.surface_areas)

    @property
    def feature_total_perimeter(self):
        return np.sum(self.perimeter)

    @property
    def feature_mean_item_perimeter(self):
        return np.mean(self.perimeter)

    @property
    def feature_field_area(self):  # todo: not defined for small n
        return self._ch.convex_hull.volume

    @property
    def feature_numerosity(self):
        return len(self._xy)

    @property
    def feature_converage(self):
        """ percent coverage in the field area. It takes thus the item size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity

        """

        try:
            return self.feature_total_surface_area / self.feature_field_area
        except:
            return None

    @property
    def feature_logSize(self):
        return misc.log2(self.feature_total_surface_area) + misc.log2(
            self.feature_mean_item_surface_area)

    @property
    def feature_logSpacing(self):
        return misc.log2(self.feature_field_area) + misc.log2(
            self.feature_sparsity)

    @property
    def feature_sparsity(self):
        return self.feature_field_area / self.feature_numerosity

    @property
    def convex_hull_positions_full(self):  # FIXME not defined for l<3
        """this convex_hull takes into account the dot diameter"""
        return self._ch.full_xy

    @property
    def feature_field_area_full(self):  # FIXME not used (correct?)
        return self._ch.full_field_area

    @property
    def object_id(self):
        """md5_hash of position, diameter"""

        m = md5()
        m.update(self._xy.tobytes())  # to byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        m.update(self.surface_areas.tobytes())
        return m.hexdigest()

    def as_dict(self, rounded_values=False):
        if rounded_values:
            xy = np.round(self._xy).astype(np.int).tolist()
            dia = np.round(self._diameters).astype(np.int).tolist()
        else:
            xy = self._xy.tolist()
            dia = self._diameters.tolist()
        return {"object_id": self.object_id, "xy": xy, "diameters": dia}

    def save(self, json_filename, rounded_values=False, indent=None):

        with open(json_filename, 'w') as fl:
            json.dump(self.as_dict(rounded_values), fl, indent=indent)

    def load(self, json_filename):

        with open(json_filename, 'r') as fl:
            dict = json.load(fl)
        self.read_from_dict(dict)

    def read_from_dict(self, dict):
        """read Dot collection from dict"""
        self._xy = np.array(dict["xy"])
        self._diameters = np.array(dict["diameters"])
        self.set_array_modified()

    def get_features_dict(self):
        """ordered dictionary with the most important feature"""
        rtn = [("Object_id", self.object_id),
               ("Numerosity", self.feature_numerosity),
               (vf.TotalSurfaceArea.label, self.feature_total_surface_area),
               (vf.ItemSurfaceArea.label, self.feature_mean_item_surface_area),
               (vf.ItemDiameter.label, self.feature_mean_item_diameter),
               (vf.ItemPerimeter.label, self.feature_mean_item_perimeter),
               (vf.TotalPerimeter.label, self.feature_total_perimeter),
               (vf.FieldArea.label, self.feature_field_area),
               (vf.Sparsity.label, self.feature_sparsity),
               (vf.Coverage.label, self.feature_converage),
               (vf.LogSize.label, self.feature_logSize),
               (vf.LogSpacing.label, self.feature_logSpacing)]
        return OrderedDict(rtn)

    def get_features_text(self, with_object_id=True, extended_format=False, spacing_char="."):
        if extended_format:
            rtn = None
            for k, v in self.get_features_dict().items():
                if rtn is None:
                    if with_object_id:
                        rtn = "- {}\n".format(v)
                    else:
                        rtn = ""
                else:
                    if rtn == "":
                        name = "- " + k
                    else:
                        name = "  " + k
                    try:
                        value = "{0:.2f}\n".format(v)  # try rounding
                    except:
                        value = "{}\n".format(v)

                    rtn += name + (spacing_char * (22 - len(name))) + (" " * (14 - len(value))) + value
        else:
            if with_object_id:
                rtn = "id: {}".format(self.object_id)
            else:
                rtn = ""
            rtn += "n: {}, TSA: {}, ISA: {}, FA: {}, SPAR: {:.3f}, logSIZE: {:.2f}, logSPACE: {:.2f} COV: {:.2f}".format(
                self.feature_numerosity,
                int(self.feature_total_surface_area),
                int(self.feature_mean_item_surface_area),
                int(self.feature_field_area),
                self.feature_sparsity,
                self.feature_logSize,
                self.feature_logSpacing,
                self.feature_converage)
        return rtn

    def set_array_modified(self):
        self._ch = misc.EfficientConvexHullDots(self._xy, self._diameters)

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
        self.set_array_modified()

    def clear(self):
        self._xy = np.array([[]])
        self._diameters = np.array([])
        self.set_array_modified()

    def delete(self, index):
        self._xy = np.delete(self._xy, index, axis=0)
        self._diameters = np.delete(self._diameters, index)
        self.set_array_modified()

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """
        if indices is None:
            indices = list(range(self.feature_numerosity))

        return DotCollection(xy=self._xy[indices, :].copy(),
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

    @property
    def feature_mean_item_diameter(self):
        return np.mean(self._diameters)

    def get_distance_matrix(self, between_positions=False):
        """between position ignores the dot size"""
        dist = spatial.distance.cdist(self._xy, self._xy)  # matrix with all distance between all points
        if not between_positions:
            # subtract dot diameter
            radii_mtx = np.ones((self.feature_numerosity, 1)) + self._diameters[:, np.newaxis].T / 2
            dist -= radii_mtx  # for each row
            dist -= radii_mtx.T  # each each column
        return dist

    @property
    def expension(self):
        """ maximal distance between to points plus diameter of the two points"""

        dist = self.get_distance_matrix(between_positions=True)
        # add dot diameter
        radii_mtx = np.ones((self.feature_numerosity, 1)) + self._diameters[:, np.newaxis].T / 2
        dist += radii_mtx  # add to each row
        dist += radii_mtx.T  # add two each column
        return np.max(dist)

    def center_array(self):
        self._xy -= self.center_of_outer_positions
        self.set_array_modified()



class DotArray(DotCollection):

    def __init__(self, target_array_radius=None,
                 minimum_gap=2,
                 xy=None,
                 diameters=None,
                 features=None,
                 dot_array_file=None):
        """Dot array is restricted to a certain area and can shuffle positions
        and find random free space and be realigned itself

        target_array_radius or dot_array_file needs to be define."""

        DotCollection.__init__(self)
        if dot_array_file is None:
            if target_array_radius is None:
                raise RuntimeError("target_array_radius need to be defined, "
                                   "if DotArray is not loaded from file.")
            self._attributes = ItemAttributesList()
            self.target_array_radius = target_array_radius
            self.minimum_gap = minimum_gap
            if xy is not None or diameters is not None or features is not None:
                self.append(xy, diameters, attributes=features)
        else:
            self.load(json_filename=dot_array_file)

    @property
    def attributes(self):
        return self._attributes

    def clear(self):
        DotCollection.clear(self)
        self._attributes.clear()

    def append(self, xy, item_diameters, attributes=None):
        """append dots using numpy array
        attributes ItemAttributes of ItemAttributesList"""

        item_diameters = misc.numpy_vector(item_diameters)
        super().append(xy=xy, item_diameters=item_diameters)

        if attributes is None:
            attributes = [ItemAttributes()] * len(item_diameters)
        elif len(item_diameters) > 1 and isinstance(attributes, ItemAttributes):
            attributes = [attributes] * len(item_diameters)

        self._attributes.append(attributes=attributes)

        if (self._attributes.length != len(self._diameters)):
            raise RuntimeError(u"Bad shaped data: " + u"attributes have not "
                                                      u"the same length as diameter")

    def delete(self, index):
        DotCollection.delete(self, index)
        self._attributes.delete(index)

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """
        if indices is None:
            indices = list(range(self.feature_numerosity))

        return DotArray(target_array_radius=self.target_array_radius,
                        minimum_gap=self.minimum_gap,
                        xy=self._xy[indices, :].copy(),
                        diameters=self._diameters[indices].copy(),
                        features=self._attributes.copy())

    def append_dot(self, dot):
        self.append(xy=[dot.x, dot.y], item_diameters=dot.diameter, attributes=dot.attributes)

    def join(self, dot_array, realign=False):
        """add another dot arrays"""

        self.append(xy=dot_array._xy, item_diameters=dot_array._diameters, attributes=dot_array.attributes)
        if realign:
            self.realign()

    def get_dots(self, indices=None, diameter=None, colour=None, picture=None):  # todo: search by attributes
        """returns all dots
         filtering possible:
         if diameter/colour/picture is defined it returns only dots a particular diameter/colour/picture

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
                    (colour is not None and att.colour != colour) or \
                    (picture is not None and att.picture != picture):
                continue

            rtn.append(Dot(x=xy[0], y=xy[1], diameter=dia, attributes=att))
        return rtn

    def _jitter_identical_positions(self, jitter_size=0.1):
        """jitters points with identical position"""

        for idx, ref_dot in enumerate(self._xy):
            identical = np.where(np.all(np.equal(self._xy, ref_dot), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != idx:
                        self._xy[x, :] -= misc.polar2cartesian([[jitter_size,
                                                                 random.random() * 2 * np.pi]])[0]

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

            tmp_polar = misc.cartesian2polar(self._xy[idx, :] - self._xy[dot_id, :])
            tmp_polar[:, 0] = 0.000000001 + minimum_gap - dist[idx]  # determine movement size
            xy = misc.polar2cartesian(tmp_polar)
            self._xy[idx, :] = np.array([self._xy[idx, 0] + xy[:, 0], self._xy[idx, 1] + xy[:, 1]]).T
            shift_required = True

        return shift_required

    def remove_overlap_from_inner_to_outer(self, minimum_gap):

        shift_required = False
        # from inner to outer remove overlaps
        for i in np.argsort(misc.cartesian2polar(self._xy, radii_only=True)):
            if self._remove_overlap_for_dot(dot_id=i, minimum_gap=minimum_gap):
                shift_required = True

        if shift_required:
            self.set_array_modified()

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

        # sqeeze in points that pop out of the stimulus area radius
        cnt = 0
        while True:
            radii = misc.cartesian2polar(self._xy, radii_only=True)
            too_far = np.where((radii + self._diameters // 2) > self.target_array_radius)[0]  # find outlier
            if len(too_far) > 0:

                # squeeze in outlier
                polar = misc.cartesian2polar([self._xy[too_far[0], :]])[0]
                polar[0] = self.target_array_radius - self._diameters[
                    too_far[0]] // 2 - 0.000000001  # new radius #todo check if 0.00001 required
                new_xy = misc.polar2cartesian([polar])[0]
                self._xy[too_far[0], :] = new_xy

                # remove overlaps centered around new outlier position
                self._xy -= new_xy
                # remove all overlaps (inner to outer, i.e. starting with outlier)
                self.remove_overlap_from_inner_to_outer(minimum_gap=self.minimum_gap)
                # new pos for outlier
                self._xy += new_xy  # move back to old position
                shift_required = True
            else:
                break  # end while loop

            cnt += 1
            if cnt > 20:
                error = True
                break

        if error:
            return False, u"Can't find solution when removing outlier (n=" + \
                   str(self.feature_numerosity) + ")"

        self.set_array_modified()
        if not shift_required:
            return True, ""
        else:
            return self.realign()  # recursion

    @property
    def feature_coverage_target_area(self):
        """density takes into account the full possible target area (i.e., stimulus radius) """
        return np.pi * self.target_array_radius ** 2 / self.feature_total_surface_area

    def split_array_by_colour(self):
        """returns a list of arrays
        each array contains all dots of with particular colour"""
        rtn = []
        for c in np.unique(self._attributes.colours):
            if c is not None:
                da = DotArray(target_array_radius=self.target_array_radius,
                              minimum_gap=self.minimum_gap)
                for d in self.get_dots(colour=c):
                    da.append_dot(d)
                rtn.append(da)
        return rtn

    def get_features_split_by_colours(self):
        """returns is unicolor or no color"""
        if len(np.unique(self._attributes.colours)) == 1:
            return None

        dicts = []
        for da in self.split_array_by_colour():
            feat = da.get_features_dict()
            feat["Object_id"] += str(da._attributes.colours[0])
            dicts.append(feat)
        return misc.join_dict_list(dicts)

    def get_csv(self, num_format="%7.2f", variable_names=True,
                object_id_column=True, num_idx_column=True,
                colour_column=False, picture_column=False):  # todo print features
        """Return the dot array as csv text

        Parameter
        ---------
        variable_names : bool, optional
            if True variable name will be printed in the first line

        """

        rtn = ""
        if variable_names:
            if object_id_column:
                rtn += u"object_id,"
            if num_idx_column:
                rtn += u"num_id,"
            rtn += u"x,y,diameter"
            if colour_column:
                rtn += u",colour"
            if picture_column:
                rtn += u",file"
            rtn += u"\n"

        obj_id = self.object_id
        for cnt in range(len(self._xy)):
            if object_id_column:
                rtn += "{0}, ".format(obj_id)
            if num_idx_column:
                rtn += "{},".format(self.feature_numerosity)
            rtn += num_format % self._xy[cnt, 0] + "," + num_format % self._xy[cnt, 1] + "," + \
                   num_format % self._diameters[cnt]
            if colour_column:
                rtn += ", {}".format(self._attributes.colours[cnt])
            if picture_column:
                rtn += ", {}".format(self._attributes.pictures[cnt])
            rtn += "\n"
        return rtn

    def random_free_dot_position(self, dot_diameter,
                                 allow_overlapping=False,
                                 prefer_inside_field_area=False,
                                 squared_array = False,
                                 occupied_space=None): #TODO rounded values
        """returns a available random xy position

        raise exception if not found
        occupied space: see generator generate
        """

        try_out_inside_convex_hull = 1000
        if prefer_inside_field_area:
            delaunay = spatial.Delaunay(self.convex_hull_positions)
        else:
            delaunay = None
        cnt = 0

        target_radius = self.target_array_radius - (dot_diameter / 2.0)
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

    def shuffle_all_positions(self, allow_overlapping=False):
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        new_diameters = np.array([])
        new_xy = np.array([])

        for d in self._diameters:
            try:
                xy = self.random_free_dot_position(d, allow_overlapping=allow_overlapping)
            except:
                raise RuntimeError("Can't shuffle dot array. No free positions.")
            new_diameters = np.append(new_diameters, d)
            if len(new_xy) == 0:
                new_xy = np.array([xy])
            else:
                new_xy = np.append(new_xy, [xy], axis=0)

        self._diameters = new_diameters
        self._xy = new_xy
        self.set_array_modified()

    def number_deviant(self, change_numerosity, prefer_keeping_field_area=False):
        """number deviant
        """

        try_out = 100
        # make a copy for the deviant
        deviant = self.copy()
        if self.feature_numerosity + change_numerosity <= 0:
            deviant.clear()
        else:
            # add or remove random dots
            for _ in range(abs(change_numerosity)):
                if prefer_keeping_field_area:
                    ch = deviant.convex_hull_indices
                else:
                    ch = []
                for x in range(try_out):
                    rnd = random.randint(0, deviant.feature_numerosity-1) # do not use np.random
                    if rnd not in ch or change_numerosity > 0:
                        break

                if change_numerosity < 0:
                    # remove dots
                    deviant.delete(rnd)
                else:
                    # copy a random dot
                    try:
                        deviant.append(xy=deviant.random_free_dot_position(dot_diameter=deviant._diameters[rnd],
                                                                           prefer_inside_field_area=prefer_keeping_field_area),
                                       item_diameters=deviant._diameters[rnd],
                                       attributes=deviant._attributes[rnd])
                    except:
                        # no free position
                        raise RuntimeError("Can't make the deviant. No free position")
        return deviant

    def as_dict(self, rounded_values=False):
        d = super().as_dict(rounded_values)
        d.update({"minimum_gap": self.minimum_gap,
             "target_array_radius": self.target_array_radius,
             "attributes": self._attributes.as_dict()})
        return d

    def read_from_dict(self, dict):
        super().read_from_dict(dict)
        self.minimum_gap = dict["minimum_gap"]
        self.target_array_radius = dict["target_array_radius"]
        self._attributes = ItemAttributesList()
        self._attributes.read_from_dict(dict["attributes"])

    def match(self, match_feature, match_dot_array=None):
        """
        match_properties: continuous property or list of continuous properties
        several properties to be matched

        if match dot array is specified, array will be match to match_dot_array, otherwise
        the values defined in match_features is used.

        some matching requires realignement to avoid overlaps. However,
        realigment might result in a different field area. Thus, realign after
        matching for  Size parameter and realign before matching space
        parameter.

        """

        # type check
        assert isinstance(match_feature, vf.ALL_VISUAL_FEATURES)

        vf.check_feature_list([match_feature],
                              check_set_value=match_dot_array
                                                             is None)

        # copy and change values to match this stimulus
        feat = copy(match_feature)
        if match_dot_array is not None:
            feat.adapt_value(match_dot_array)

        # Adapt
        if isinstance(feat, vf.ItemDiameter):
            self._match_item_diameter(mean_item_diameter=feat.value)

        elif isinstance(feat, vf.ItemPerimeter):
            self._match_item_diameter(mean_item_diameter=feat.value/np.pi)

        elif isinstance(feat, vf.TotalPerimeter):
            mean_dot_diameter = feat.value / (self.feature_numerosity * np.pi)
            self._match_item_diameter(mean_dot_diameter)

        elif isinstance(feat, vf.ItemSurfaceArea):
            ta = self.feature_numerosity * feat.value
            self._match_total_surface_area(surface_area=ta)

        elif isinstance(feat, vf.TotalSurfaceArea):
            self._match_total_surface_area(surface_area=feat.value)

        elif isinstance(feat, vf.LogSize):
            logtsa = 0.5 * feat.value + 0.5 * misc.log2(self.feature_numerosity)
            self._match_total_surface_area(2 ** logtsa)

        elif isinstance(feat, vf.LogSpacing):
            logfa = 0.5 * feat.value + 0.5 * misc.log2(self.feature_numerosity)
            self._match_field_area(field_area=2 ** logfa,
                                   precision=feat.spacing_precision)

        elif isinstance(feat, vf.Sparsity):
            fa = feat.value * self.feature_numerosity
            self._match_field_area(field_area=fa,
                                   precision=feat.spacing_precision)

        elif isinstance(feat, vf.FieldArea):
            self._match_field_area(field_area=feat.value,
                                   precision=feat.spacing_precision)

        elif isinstance(feat, vf.Coverage):
            self._match_coverage(coverage=feat.value,
                                 precision=feat.spacing_precision,
                                 match_FA2TA_ratio=feat.match_ratio_fieldarea2totalarea) #FIXME experimemtal


    def _match_total_surface_area(self, surface_area):
        # changes diameter
        a_scale = (surface_area / self.feature_total_surface_area)
        self._diameters = np.sqrt(self.surface_areas * a_scale) * 2 / np.sqrt(
            np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
        self.set_array_modified()

    def _match_item_diameter(self, mean_item_diameter):
        # changes diameter

        scale = mean_item_diameter / self.feature_mean_item_diameter
        self._diameters = self._diameters * scale
        self.set_array_modified()

    def _match_field_area(self, field_area,
                          precision=vf._DEFAULT_SPACING_PRECISION,
                          use_scaling_only=False):
        """changes the convex hull area to a desired size with certain precision

        uses scaling radial positions if field area has to be increased
        uses replacement of outer points (and later re-scaling)

        iterative method can takes some time.
        """

        if self.feature_field_area is None:
            return  # not defined
        elif field_area > self.feature_field_area or use_scaling_only:
            # field area is currently too small or scaling is enforced
            return self.__scale_field_area(field_area=field_area,
                                           precision=precision)
        elif field_area < self.feature_field_area:
            # field area is too large
            self.__decrease_field_area_by_replacement(max_field_area=field_area)
            # ..and rescaling to avoid to compensate for possible too
            # strong decrease
            return self.__scale_field_area(field_area=field_area,
                                           precision=precision)
        else:
            return

    def __decrease_field_area_by_replacement(self, max_field_area):
        """decreases filed area by recursively moving the most outer point
        to some more central free position (avoids overlapping)

        return False if not possible else True"""


        # centered points
        old_center = self.center_of_outer_positions
        self._xy -= old_center

        removed_dots = []
        while self.feature_field_area > max_field_area:
            # remove one random outer dot and remember it
            vertices = self._ch.convex_hull.vertices
            idx = vertices[random.randint(0, len(vertices)-1)]

            removed_dots.extend(self.get_dots(indices=[idx]))
            self.delete(idx)

        # add dots to free pos inside the convex hall
        for d in removed_dots:
            d.xy = self.random_free_dot_position(d.diameter,
                                               allow_overlapping=False,
                                               prefer_inside_field_area=True)
            self.append_dot(d)

        self._xy += old_center
        self.set_array_modified()

    def __scale_field_area(self, field_area,
                           precision=vf._DEFAULT_SPACING_PRECISION):
        """change the convex hull area to a desired size by scale the polar
        positions  with certain precision

        iterative method can takes some time.
        """
        current = self.feature_field_area

        if current is None:
            return  # not defined

        scale = 1  # find good scale
        step = 0.1
        if field_area < current:  # current too larger
            step *= -1

        # centered points
        old_center = self.center_of_outer_positions
        self._xy -= old_center
        centered_polar = misc.cartesian2polar(self._xy)

        # iteratively determine scale
        while abs(current - field_area) > precision:

            scale += step

            self._xy = misc.polar2cartesian(centered_polar * [scale, 1])
            self.set_array_modified()  # required to recalc convex hull
            current = self.feature_field_area

            if (current < field_area and step < 0) or \
                    (current > field_area and step > 0):
                step *= -0.2  # change direction and finer grain

        self._xy += old_center
        self.set_array_modified()

    def _match_coverage(self, coverage, precision=vf._DEFAULT_SPACING_PRECISION,
                        match_FA2TA_ratio=0.5):
        # FIXME check drifting outwards if extra space is small and match_FA2TA_ratio=1
        # FIXME when to realign, realignment changes field_area!
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """

        print("WARNING: _match_coverage is a experimental ")
        # dens = convex_hull_area / total_surface_area
        if match_FA2TA_ratio < 0 or match_FA2TA_ratio > 1:
            match_FA2TA_ratio = 0.5

        total_area_change100 = (coverage * self.feature_field_area) - self.feature_total_surface_area
        d_change_total_area = total_area_change100 * (1 - match_FA2TA_ratio)
        if abs(d_change_total_area) > 0:
            self._match_total_surface_area(surface_area=self.feature_total_surface_area + d_change_total_area)

        self._match_field_area(field_area=self.feature_total_surface_area / coverage,
                               precision=precision)
