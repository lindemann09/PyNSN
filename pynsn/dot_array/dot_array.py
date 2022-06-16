"""
Dot Array
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import random
import numpy as np
from scipy import spatial

from .cloud import _Cloud
from ..lib import misc, geometry
from pynsn.dot_array.shape import Dot

# TODO: How to deal with rounding? Is saving to precises? Suggestion:
#  introduction precision parameter that is used by as_dict and get_csv and
#  hash

class DotArray(_Cloud):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    def __init__(self, target_array_radius, minimum_gap,
                 xy=None, diameters=None, attributes=None):
        """Dot array is restricted to a certain area, it has a target area
        and a minimum gap.
        This features allow find shuffling free position and matching
        features.

        target_array_radius or dot_array_file needs to be define."""

        super().__init__(xy=xy, attributes=attributes)
        if diameters is None:
            self._diameters = np.array([])
        else:
            self._diameters = misc.numpy_vector(diameters)
        if self._xy.shape[0] != len(self._diameters):
            raise RuntimeError("Bad shaped data: " +
                    u"xy has not the same length as item_diameters")

        self.target_array_radius = target_array_radius
        self.minimum_gap = minimum_gap

    @property
    def diameters(self):
        return self._diameters

    @property
    def surface_areas(self):
        # a = pi r**2 = pi d**2 / 4
        return np.pi * (self._diameters ** 2) / 4.0

    @property
    def perimeter(self):
        return np.pi * self._diameters

    def round(self, decimals=0, int_type=np.int64):
        """Round values of the array."""

        if decimals is None:
            return
        super().round(decimals, int_type)
        self._diameters = np.round(self._diameters, decimals=decimals)
        if decimals == 0:
            self._diameters = self._diameters.astype(int_type)

    def as_dict(self):
        """
        """
        d = super().as_dict()
        d.update({"diameters": self._diameters.tolist(),
                  "minimum_gap": self.minimum_gap,
                  "target_array_radius": self.target_array_radius})
        return d

    def read_from_dict(self, dict):
        """read Dot collection from dict"""
        super().read_from_dict(dict)
        self._diameters = np.array(dict["diameters"])
        self.minimum_gap = dict["minimum_gap"]
        self.target_array_radius = dict["target_array_radius"]

    def add(self, dots):
        """append one dot or list of dots"""
        if not isinstance(dots, (list, tuple)):
            dots = [dots]
        for d in dots:
            assert isinstance(d, Dot)
            self._append_xy_attribute(xy=d.xy, attributes=d.attribute)
            self._diameters = np.append(self._diameters, d.diameter)

    def clear(self):
        super().clear()
        self._diameters = np.array([])

    def delete(self, index):
        super().delete(index)
        self._diameters = np.delete(self._diameters, index)

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if indices is None:
            indices = list(range(self._features.numerosity))

        return DotArray(target_array_radius=self.target_array_radius,
                        minimum_gap=self.minimum_gap,
                        xy=self._xy[indices, :].copy(),
                        diameters=self._diameters[indices].copy(),
                        attributes=self._attributes[indices].copy())

    def distances(self, xy, diameter):
        """Distances toward a single point (numpy arrays xy, diameter)

        negative numbers indicate overlap
        """
        if len(self._xy) == 0:
            return np.array([])
        else:
            return np.hypot(self._xy[:, 0] - xy[0], self._xy[:, 1] - xy[1]) - \
                   ((self._diameters + diameter) / 2.0)

    @property
    def center_of_mass(self):
        weighted_sum = np.sum(self._xy * self._diameters[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self._diameters)

    def get(self, indices=None):
        """returns all dots

        indices int or list of ints
         filtering possible:
         if diameter is defined it returns only dots a particular
         diameter

         """

        if indices is None:
            return [Dot(x=xy[0], y=xy[1], diameter=dia, attribute=att) \
                    for xy, dia, att in zip(self._xy, self._diameters,
                                            self._attributes)]
        try:
            indices = list(indices)  # check if iterable
        except:
            indices = [indices]

        return [Dot(x=xy[0], y=xy[1], diameter=dia, attribute=att) \
                    for xy, dia, att in zip(self._xy[indices, :],
                                        self._diameters[indices],
                                        self._attributes[indices])]

    def find(self, diameter=None, item_attributes=None):
        """returns all dots
        indices : int or list of ins

         filtering possible:
         if diameter/colour/picture is defined it returns only dots a
         particular diameter/colour/picture

         """
        rtn = []
        for xy, dia, att in zip(self._xy, self._diameters, self._attributes):
            if (diameter is not None and dia != diameter) or \
                    (item_attributes is not None and att !=  item_attributes):
                continue

            rtn.append(Dot(x=xy[0], y=xy[1], diameter=dia, attribute=att))
        return rtn

    def csv(self, variable_names=True,
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
                rtn += "{},".format(self._features.numerosity)
            #rtn += num_format % self._xy[cnt, 0] + "," + num_format %
                # self._xy[cnt, 1] + "," + \
            #       num_format % self._diameters[cnt]
            rtn += "{},{},{}".format(self._xy[cnt, 0], self._xy[cnt, 1],
                                     self._diameters[cnt])
            if attribute_column:
                rtn += ", {}".format(self._attributes[cnt])
            rtn += "\n"
        return rtn

    def join(self, dot_array):
        """add another dot arrays"""
        assert isinstance(dot_array, DotArray)
        self.add(dot_array.get())

    def _remove_overlap_from_inner_to_outer(self, minimum_gap):

        shift_required = False
        # from inner to outer remove overlaps
        for i in np.argsort(geometry.cartesian2polar(self._xy, radii_only=True)):
            dist = self.distances(self._xy[i, :], self._diameters[i])
            idx_overlaps = np.where(dist < minimum_gap)[0].tolist()  # overlapping dot ids
            if len(idx_overlaps) > 1:
                shift_required = True
                idx_overlaps.remove(i)  # don't move yourself
                replace_size = minimum_gap - dist[idx_overlaps] # dist is mostly negative, because of overlap
                self._radial_replacement_from_reference_dots(ref_pos_id=i,
                                            neighbour_ids=idx_overlaps,
                                            replacement_size= replace_size)

        if shift_required:
            self._features.reset()
        return shift_required

    def realign(self): #TODO update realignment
        """Realigns the dots in order to remove all dots overlaps and dots
        outside the target area.

        If two dots overlap, the dots that is further apart from the array
        center will be moved opposite to the direction of the other dot until
        there is no overlap (note: minimum_gap parameter). If two dots have
        exactly the same position the same position one is minimally shifted
        in a random direction.

        Note: Realigning might change the field area! Match Space parameter after
        realignment.

        """

        error = False

        shift_required = self._remove_overlap_from_inner_to_outer(minimum_gap=self.minimum_gap)

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
                self._remove_overlap_from_inner_to_outer(minimum_gap=self.minimum_gap)
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
                   str(self._features.numerosity) + ")"

        self._features.reset()
        if not shift_required:
            return True, ""
        else:
            return self.realign()  # recursion

    def random_free_position(self, dot_diameter,
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
            delaunay = spatial.Delaunay(self._features.convex_hull._xy)
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
            #Note! np.random generates identical numbers under multiprocessing

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
                xy = self.random_free_position(d,
                                               allow_overlapping=allow_overlapping,
                                               min_distance_area_boarder=min_distance_area_boarder)
            except:
                raise RuntimeError("Can't shuffle dot array. No free positions.")

            if new_xy is None:
                new_xy = np.array([xy])
            else:
                new_xy = np.append(new_xy, [xy], axis=0)

        self._xy = new_xy
        self._features.reset()

    def number_deviant(self, change_numerosity, prefer_keeping_field_area=False):
        """number deviant
        """

        TRY_OUT = 100
        # make a copy for the deviant
        deviant = self.copy()
        if self._features.numerosity + change_numerosity <= 0:
            deviant.clear()
        else:
            # add or remove random dots
            rnd =None
            for _ in range(abs(change_numerosity)):
                if prefer_keeping_field_area:
                    ch = deviant._features.convex_hull.indices
                else:
                    ch = []
                for x in range(TRY_OUT):##
                    rnd = random.randint(0, deviant._features.numerosity - 1) # do not use np.random
                    if rnd not in ch or change_numerosity > 0:
                        break

                if change_numerosity < 0:
                    # remove dots
                    deviant.delete(rnd)
                else:
                    # copy a random dot
                    rnd_dot = self.get([rnd])[0]
                    try:
                        rnd_dot.xy = deviant.random_free_position(
                            dot_diameter=rnd_dot.diameter,
                            prefer_inside_field_area=prefer_keeping_field_area)
                    except:
                        # no free position
                        raise RuntimeError("Can't make the deviant. No free position")
                    deviant.add(rnd_dot)

        return deviant

    def split_array_by_attributes(self):
        """returns a list of arrays
        each array contains all dots of with particular colour"""
        att = self._attributes
        att[np.where(att==None)] = "None" # TODO check "is none"

        rtn = []
        for c in np.unique(att):
            if c is not None:
                da = DotArray(target_array_radius=self.target_array_radius,
                              minimum_gap=self.minimum_gap)
                da.add(self.find(item_attributes =c))
                rtn.append(da)
        return rtn