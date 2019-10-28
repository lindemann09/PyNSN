"""
Dot Array
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from copy import copy
import random
from multiprocessing import Pool

import numpy as np
from scipy import spatial

from .geometry import Dot
from . import misc
from .dot_collection import DotCollection
from .item_attributes import ItemAttributesList, ItemAttributes
from . import visual_features

TWO_PI = 2 * np.pi


class DotArray(DotCollection):

    def __init__(self, target_array_radius=None,
                 minimum_gap=1,
                 xy=None,
                 diameters=None,
                 features=None,
                 dot_array_file=None):
        """Dot array is restricted to a certain area and can generate random dots and
            be realigned

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
                        self._xy[x, :] -= misc.polar2cartesian([[jitter_size, random.random() * TWO_PI]])[0]

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

    def realign(self, center_array=False):
        """Realigns the dots in order to remove all dots overlaps. If two dots
        overlap, the dots that is further apart from the arry center will be
        moved opposite to the direction of the other dot until there is no
        overlap (note: minimun_gap parameter). If two dots have exactly the same
        position the same position one is minimally shifted in a random direction.

        """

        error = False

        if center_array:
            self._xy -= self.center_of_outer_positions
            self.set_array_modified()

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
                # new pos for outlyer
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
        if not shift_required:
            return True, ""
        else:
            self.set_array_modified()
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

    def __str__(self):
        return self.get_csv()

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
                                 ignore_overlapping=False,
                                 prefer_inside_field_area=False,
                                 squared_array = False,
                                 occupied_space=None):
        """moves a dot to an available random position

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
            #  proposal_polar =  np.random.random(2) * (target_radius, TWO_PI)
            #proposal_xy = misc.polar2cartesian([proposal_polar])[0]
            proposal_xy = np.random.random(2) * 2 * self.target_array_radius - \
                        self.target_array_radius
            bad_position = False
            if not squared_array:
                bad_position =  target_radius <= \
                                np.hypot(proposal_xy[0], proposal_xy[1])

            if not bad_position and prefer_inside_field_area and cnt < \
                    try_out_inside_convex_hull:
                bad_position = delaunay.find_simplex(proposal_xy) < 0

            if not bad_position and not ignore_overlapping:
                # find bad_positions
                dist = self.distances(proposal_xy, dot_diameter)
                if occupied_space:
                    dist = np.append(dist, occupied_space.distances(proposal_xy, dot_diameter))
                idx = np.where(dist < self.minimum_gap)[0]  # overlapping dot ids
                bad_position = len(idx) > 0

            if not bad_position:
                return proposal_xy
            elif cnt > 3000:
                raise RuntimeError(u"Can't find a free position")

    def shuffle_all_positions(self, ignore_overlapping=False):
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        new_diameters = np.array([])
        new_xy = np.array([])

        for d in self._diameters:
            try:
                xy = self.random_free_dot_position(d, ignore_overlapping=ignore_overlapping)
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
                    rnd = random.randint(0,
                                         deviant.feature_numerosity - 1)  # strange: do not use np.rand, because parallel running process produce the same numbers
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

    def match(self, match_features, center_array=True,
              realign = False,
              match_dot_array=None):
        """
        match_properties: continuous property or list of continuous properties
        several properties to be matched

        if match dot array is specified, array will be match to match_dot_array, otherwise
        the values defined in match_features is used.
        """

        # type check
        if not isinstance(match_features, (list, tuple)):
            match_features = [match_features]

        visual_features.check_feature_list(match_features, check_set_value=match_dot_array is None)

        # copy and change values to match this stimulus
        if match_dot_array is None:
            match_feat = match_features
        else:
            match_feat = []
            for m in match_features:
                m = copy(m)
                m.adapt_value(match_dot_array)
                match_feat.append(m)

        # Adapt
        for feat in match_feat:
            if isinstance(feat, visual_features.ItemDiameter):
                self._match_item_diameter(mean_item_diameter=feat.value)

            elif isinstance(feat, visual_features.ItemPerimeter):
                self._match_item_diameter(mean_item_diameter=feat.value/np.pi)

            elif isinstance(feat, visual_features.TotalPerimeter):
                mean_dot_diameter = feat.value / (self.feature_numerosity * np.pi)
                self._match_item_diameter(mean_dot_diameter)

            elif isinstance(feat, visual_features.ItemSurfaceArea):
                ta = self.feature_numerosity * feat.value
                self._match_total_surface_area(surface_area=ta)

            elif isinstance(feat, visual_features.TotalSurfaceArea):
                self._match_total_surface_area(surface_area=feat.value)

            elif isinstance(feat, visual_features.LogSize):
                logtsa = 0.5 * feat.value + 0.5 * misc.log2(self.feature_numerosity)
                self._match_total_surface_area(2 ** logtsa)

            elif isinstance(feat, visual_features.LogSpacing):
                logfa = 0.5 * feat.value + 0.5 * misc.log2(self.feature_numerosity)
                self._match_field_area(field_area=2 ** logfa,
                                       precision=feat.spacing_precision)

            elif isinstance(feat, visual_features.Sparsity):
                fa = feat.value * self.feature_numerosity
                self._match_field_area(field_area=fa,
                                       precision=feat.spacing_precision)

            elif isinstance(feat, visual_features.FieldArea):
                self._match_field_area(field_area=feat.value,
                                       precision=feat.spacing_precision)

            elif isinstance(feat, visual_features.Coverage):
                self._match_coverage(coverage=feat.value,
                                     precision=feat.spacing_precision,
                                     match_FA2TA_ratio=feat.match_ratio_fieldarea2totalarea)

        if realign:
            self.realign(center_array=center_array)
        elif center_array:
            self._xy -= self.center_of_outer_positions
            self.set_array_modified()

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

    def _match_field_area(self, field_area, precision=visual_features._DEFAULT_SPACING_PRECISION):
        """changes the convex hull area to a desired size with certain precision

        iterative method can takes some time.
        """
        current = self.feature_field_area

        if current is None:
            return  # not defined

        # iteratively determine scale
        scale = 1  # find good scale
        step = 0.1
        if field_area < current:  # current too larger
            step *= -1

        # centered points
        old_center = self.center_of_outer_positions
        self._xy -= old_center
        centered_polar = misc.cartesian2polar(self._xy)

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

    def _match_coverage(self, coverage, precision=visual_features._DEFAULT_SPACING_PRECISION,
                        match_FA2TA_ratio=0.5):  # FIXME check drifting outwards if extra space is small and match_FA2TA_ratio=1
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """

        # dens = convex_hull_area / total_surface_area
        if match_FA2TA_ratio < 0 or match_FA2TA_ratio > 1:
            match_FA2TA_ratio = 0.5

        total_area_change100 = (coverage * self.feature_field_area) - self.feature_total_surface_area
        d_change_total_area = total_area_change100 * (1 - match_FA2TA_ratio)
        if abs(d_change_total_area) > 0:
            self._match_total_surface_area(surface_area=self.feature_total_surface_area + d_change_total_area)

        self._match_field_area(field_area=self.feature_total_surface_area / coverage,
                               precision=precision)




class DotArrayFactory(object):

    def __init__(self,
                 target_area_radius,
                 item_diameter_mean,
                 item_diameter_range=None,
                 item_diameter_std=None,
                 item_colour=None,  # todo feature
                 minimum_gap=1):  # TODO check minim gap

        """Specification of a Random Dot Array

        Parameters:
        -----------
        stimulus_area_radius : int
            the radius of the stimulus area
        n_dots : int
            number of moving dots

        automatic logging log only the create process. If colours a changes later they are not log.
        Use manual logging in this case.

        """

        if item_diameter_std <= 0:
            item_diameter_std = None
        elif item_diameter_range is not None and \
                (item_diameter_mean <= item_diameter_range[0] or
                 item_diameter_mean >= item_diameter_range[1] or
                 item_diameter_range[0] >= item_diameter_range[1]):
            raise RuntimeError("item_diameter_mean has to be inside the defined item_diameter_range")

        self.minimum_gap = minimum_gap
        self.target_array_radius = target_area_radius
        self.item_diameter_range = item_diameter_range
        self.item_diameter_mean = item_diameter_mean
        self.item_diameter_std = item_diameter_std
        self.item_attributes = ItemAttributes(colour=item_colour)

    def generate(self, n_dots, occupied_space=None, logger=None):
        """occupied_space is a dot array (used for multicolour dot array (join after)

        returns None if not possible
        """

        rtn = DotArray(target_array_radius=self.target_array_radius,  # - distance_field_edge ?
                       minimum_gap=self.minimum_gap)

        for _ in range(n_dots):

            # diameter
            if self.item_diameter_range is None or self.item_diameter_std is None:
                # constant diameter
                diameter = self.item_diameter_mean
            else:
                # draw diameter from beta distribution
                parameter = misc.shape_parameter_beta(self.item_diameter_range,
                                                      self.item_diameter_mean,
                                                      self.item_diameter_std)
                diameter = misc.random_beta(
                    self.item_diameter_range, parameter)

            try:
                xy = rtn.random_free_dot_position(dot_diameter=diameter, occupied_space=occupied_space)
            except:
                return None
            rtn.append(xy=xy, item_diameters=diameter, attributes=self.item_attributes)

        if logger is not None:
            from .logging import LogFile # to avoid circular import
            if not isinstance(logger, LogFile): #
                raise RuntimeError("logger has to be None or a GeneratorLogger")
            logger.log(rtn)

        return rtn

    def as_dict(self):
        return {"target_array_radius": self.target_array_radius,
                "dot_diameter_mean": self.item_diameter_mean,
                "dot_diameter_range": self.item_diameter_range,
                "dot_diameter_std": self.item_diameter_std,
                "dot_colour": self.item_attributes.colour.colour,  ##todo feature
                "minimum_gap": self.minimum_gap}

    def generate_iter(self, list_of_n_dots, occupied_space=None, logger=None,
                      multiprocessing=False):  # TODO  never checked
        args = map(lambda x: (self, x, occupied_space, logger), list_of_n_dots)

        if multiprocessing:
            return Pool().imap(_make_imap_helper, args)
        else:
            return map(_make_imap_helper, args)


def _make_imap_helper(args):
    generator = args[0]
    return generator.plot(reference_dot_array=args[1],
                          occupied_space=args[2],
                          logger=args[3])
