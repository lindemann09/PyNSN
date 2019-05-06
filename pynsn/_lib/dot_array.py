"""
Dot Array
"""
from __future__ import absolute_import, print_function, division, unicode_literals
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import random
import numpy as np
from scipy import spatial
from multiprocessing import Pool

from .dot import Dot
from . import misc
from .dot_collection import DotCollection
from .item_attributes import ItemAttributesList
from .item_attributes import ItemAttributes

TWO_PI = 2 * np.pi


class DotArray(DotCollection):

    def __init__(self, target_array_radius, minimum_gap=1, xy=None,
                 diameters=None, features=None):
        """Dot array is restricted to a certain area and can generate random dots and
            be realigned """

        DotCollection.__init__(self)
        self._attributes = ItemAttributesList()
        self.target_array_radius = target_array_radius
        self.minimum_gap = minimum_gap
        if xy is not None or diameters is not None or features is not None:
            self.append(xy, diameters, attributes=features)

    @property
    def attributes(self):
        return self._attributes

    def clear(self):
        DotCollection.clear(self)
        self._attributes.clear()

    def append(self, xy, item_diameters, attributes=None):
        """append dots using numpy array
        features None, ItemFeatureList of ItemFeatures"""

        item_diameters = misc.numpy_vector(item_diameters)
        DotCollection.append(self, xy=xy, item_diameters=item_diameters)

        if attributes is None:
            attributes = ItemAttributesList(colours=[None] * len(item_diameters))
        elif len(item_diameters) > 1 and isinstance(attributes, ItemAttributes):
            tmp = ItemAttributesList()
            tmp.append_attributes(attributes)
            attributes = tmp.repeat(len(item_diameters))

        self._attributes.append_attributes(attributes=attributes)

        if (self._attributes.length != len(self._diameters)):
            raise RuntimeError(u"Bad shaped data: " + u"colour and/or picture has not the same length as diameter")

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
        self.append(xy=[dot.x, dot.y], item_diameters=dot.diameter, attributes=dot.features)

    def join(self, dot_array, realign=False):
        """add another dot arrays"""

        self.append(xy=dot_array._xy, item_diameters=dot_array._diameters, attributes=dot_array.features)
        if realign:
            self.realign()

    def get_dots(self, indices=None, diameter=None, colour=None, picture=None):  # todo: search by features
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

        for xy, dia, feat in zip(self._xy, self._diameters, self._attributes):
            i += 1
            if (indices is not None and i not in indices) or \
                    (diameter is not None and dia != diameter) or \
                    (colour is not None and feat.colour != colour) or \
                    (picture is not None and feat.picture != picture):
                continue

            rtn.append(Dot(x=xy[0], y=xy[1], diameter=dia, features=feat))
        return rtn

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
            radii = DotCollection._cartesian2polar(self._xy, radii_only=True)
            too_far = np.where((radii + self._diameters // 2) > self.target_array_radius)[0]  # find outlier
            if len(too_far) > 0:

                # squeeze in outlier
                polar = DotCollection._cartesian2polar([self._xy[too_far[0], :]])[0]
                polar[0] = self.target_array_radius - self._diameters[
                    too_far[0]] // 2 - 0.000000001  # new radius #todo check if 0.00001 required
                new_xy = DotCollection._polar2cartesian([polar])[0]
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
            feat["object_id"] += str(da._attributes.colours[0])
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
                                 occupied_space=None):
        """moves a dot to an available random position

        raise exception if not found
        occupied space: see generator make
        """

        try_out_inside_convex_hull = 1000
        if prefer_inside_field_area:
            delaunay = spatial.Delaunay(self.convex_hull_positions)
        else:
            delaunay = None
        cnt = 0
        while True:
            cnt += 1
            proposal_polar = np.random.random(2) * \
                             (self.target_array_radius - (dot_diameter / 2.0), TWO_PI)
            proposal_xy = DotCollection._polar2cartesian([proposal_polar])[0]

            bad_position = False
            if prefer_inside_field_area and cnt < try_out_inside_convex_hull:
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

    def yaml_dump(self, document_separator=True):
        import yaml

        if document_separator:
            rtn = "---\n# {}, {}\n".format(self.object_id, self.feature_numerosity)
        else:
            rtn = ""
        d = {"minimum_gap": self.minimum_gap,
             "target_array_radius": self.target_array_radius,
             "xy": self._xy.tolist(),
             "diameter": self._diameters.tolist(),
             "features": self._attributes.as_dict()}
        return rtn + yaml.dump(d)



class RandomDotArray(object):

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
        self.item_feature = ItemAttributes(colour=item_colour)

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
            rtn.append(xy=xy, item_diameters=diameter, attributes=self.item_feature)

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
                "dot_colour": self.item_feature.colour.colour,  ##todo feature
                "minimum_gap": self.minimum_gap}

    def generate_iter(self, list_of_n_dots, occupied_space=None, logger=None, multiprocessing=False):  # TODO  never checked
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
