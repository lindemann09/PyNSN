"""
Dot Array
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np

from ._generic_object_array import GenericObjectArray
from .. import misc, geometry
from ..shapes import Dot
from . import _manipulate

# TODO: How to deal with rounding? Is saving to precises? Suggestion:
#  introduction precision parameter that is used by as_dict and get_csv and
#  hash


class DotArray(GenericObjectArray):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    def __init__(self,
                 target_area_radius,
                 min_dist_between = 2,
                 min_dist_area_boarder = 1,
                 xy = None,
                 diameters = None,
                 attributes = None):
        """Dot array is restricted to a certain area, it has a target area
        and a minimum gap.

        This features allows shuffling free position and adapting
        features.
        """
        super().__init__(xy=xy, attributes=attributes,
                         target_area_radius=target_area_radius,
                         min_dist_between=min_dist_between,
                         min_dist_area_boarder=min_dist_area_boarder)
        if diameters is None:
            self._diameters = np.array([])
        else:
            self._diameters = misc.numpy_vector(diameters)
        if self._xy.shape[0] != len(self._diameters):
            raise ValueError("Bad shaped data: " +
                             u"xy has not the same length as item_diameters")

    def add(self, dots):
        """append one dot or list of dots"""
        if not isinstance(dots, (list, tuple)):
            dots = [dots]
        for d in dots:
            assert isinstance(d, Dot)
            self._append_xy_attribute(xy=d.xy, attributes=d.attribute)
            self._diameters = np.append(self._diameters, d.diameter)

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

    def round(self, decimals=0, int_type=np.int32):
        """Round values of the array."""

        if decimals is None:
            return
        self._xy = misc.numpy_round2(self._xy, decimals=decimals,
                                     int_type=int_type)
        self._diameters = misc.numpy_round2(self._diameters, decimals=decimals,
                                            int_type=int_type)

    def as_dict(self):
        """
        """
        d = super().as_dict()
        d.update({"diameters": self._diameters.tolist(),
                  "min_dist_between": self.min_dist_between,
                  "target_area_radius": self.target_area_radius})
        return d

    def read_from_dict(self, the_dict):
        """read Dot collection from dict"""
        super().read_from_dict(the_dict)
        self._diameters = np.array(the_dict["diameters"])
        self.min_dist_between = the_dict["min_dist_between"]
        self.target_area_radius = the_dict["target_area_radius"]

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

        return DotArray(target_area_radius=self.target_area_radius,
                        min_dist_between=self.min_dist_between,
                        min_dist_area_boarder = self.min_dist_area_boarder,
                        xy=self._xy[indices, :].copy(),
                        diameters=self._diameters[indices].copy(),
                        attributes=self._attributes[indices].copy())

    def distances(self, dot):
        """Distances toward a single dot
        negative numbers indicate overlap

        Returns
        -------
        distances : numpy array of distances
        """
        assert isinstance(dot, Dot)
        if len(self._xy) == 0:
            return np.array([])
        else:
            rtn = np.hypot(self._xy[:, 0] - dot.x, self._xy[:, 1] - dot.y) - \
                   ((self._diameters + dot.diameter) / 2.0)
            return rtn

    def get(self, indices=None):
        """returns all dots

        indices int or list of ints
        """

        if indices is None:
            return [Dot(xy=xy, diameter=dia, attribute=att) \
                    for xy, dia, att in zip(self._xy, self._diameters,
                                            self._attributes)]
        try:
            indices = list(indices)  # check if iterable
        except TypeError:
            indices = [indices]

        return [Dot(xy=xy, diameter=dia, attribute=att) \
                for xy, dia, att in zip(self._xy[indices, :],
                                        self._diameters[indices],
                                        self._attributes[indices])]

    def find(self, diameter=None, attribute=None):
        """returns indices of found objects
        """
        rtn = []
        for i in range(len(self._diameters)):
            if (diameter is not None and self._diameters[i] != diameter) or \
                    (attribute is not None and self._attributes[i] != attribute):
                continue
            rtn.append(i)
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

    def get_random_free_position(self,
                                 dot_diameter,
                                 allow_overlapping = False,
                                 inside_convex_hull = False,
                                 occupied_space = None):
        """returns a available random xy position

        raise exception if not found
        occupied space: see generator generate
        """

        if isinstance(occupied_space, GenericObjectArray):
            os_distance_fnc = occupied_space.distances
        else:
            os_distance_fnc = None
        if inside_convex_hull:
            convex_hull_xy = self.features.convex_hull_positions.xy
        else:
            convex_hull_xy = None
        return _manipulate.get_random_free_position(
            the_object=Dot(xy=(0, 0), diameter=dot_diameter),
            target_area_radius=self.target_area_radius,
            allow_overlapping=allow_overlapping,
            distances_function=self.distances,
            min_dist_between=self.min_dist_between,
            min_dist_area_boarder=self.min_dist_area_boarder,
            occupied_space_distances_function=os_distance_fnc,
            convex_hull_xy=convex_hull_xy)

    def _remove_overlap_from_inner_to_outer(self):

        shift_required = False
        # from inner to outer remove overlaps
        for i in np.argsort(geometry.cartesian2polar(self._xy, radii_only=True)):
            dist_mtx = self.distances_matrix(between_positions=False,
                                             overlap_is_zero=False)

            idx_overlaps = np.where(dist < self.min_dist_between)[0].tolist()  # overlapping dot ids
            if len(idx_overlaps) > 1:
                shift_required = True
                idx_overlaps.remove(i)  # don't move yourself
                replace_size =self.min_dist_between - dist[idx_overlaps]  # dist is mostly negative, because of overlap
                self._xy = _manipulate.radial_replacement_from_reference_dots(
                                             xy=self._xy,
                                             ref_pos_id=i,
                                             neighbour_ids=idx_overlaps,
                                             replacement_size=replace_size)
        return shift_required

    def realign_old(self):  # TODO update realignment
        """Realigns the obejcts in order to remove all dots overlaps and dots
        outside the target area.

        If two dots overlap, the dots that is further apart from the array
        center will be moved opposite to the direction of the other dot until
        there is no overlap (note: min_dist_between parameter). If two dots have
        exactly the same position the same position one is minimally shifted
        in a random direction.

        Note: Realigning might change the field area! Adapt Space parameter after
        realignment.

        """

        error = False

        shift_required = self._remove_overlap_from_inner_to_outer()

        # sqeeze in points that pop out of the image area radius
        cnt = 0
        while True:
            radii = geometry.cartesian2polar(self._xy, radii_only=True)
            too_far = np.where((radii + self._diameters / 2) > self.target_area_radius)[0]  # find outlier
            if len(too_far) > 0:

                # squeeze in outlier
                polar = geometry.cartesian2polar([self._xy[too_far[0], :]])[0]
                polar[0] = self.target_area_radius - self._diameters[
                    too_far[0]] / 2 - 0.000000001  # new radius #todo check if 0.00001 required
                new_xy = geometry.polar2cartesian([polar])[0]
                self._xy[too_far[0], :] = new_xy

                # remove overlaps centered around new outlier position
                self._xy = self._xy - new_xy
                # remove all overlaps (inner to outer, i.e. starting with outlier)
                self._remove_overlap_from_inner_to_outer()
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
    
    def realign(self):
        error = False
        realign_required = False
        cnt = 0
        while True:
            cnt += 1
            if cnt > 20:
                error = True
                break

        if error:
            raise RuntimeError("Can't find solution when removing outlier (n=" + \
                   str(self._features.numerosity) + ")")

        self._features.reset()
        if not realign_required:
            return True, ""
        else:
            return self.realign()  # recursion

    def shuffle_all_positions(self, allow_overlapping=False):
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        new_xy = None
        for d in self.diameters:
            try:
                xy = self.get_random_free_position(d,
                                                   allow_overlapping=allow_overlapping)
            except StopIteration as e:
                raise StopIteration("Can't shuffle dot array. No free positions found.")

            if new_xy is None:
                new_xy = np.array([xy])
            else:
                new_xy = np.append(new_xy, [xy], axis=0)

        self._xy = new_xy
        self._features.reset()


    def get_split_arrays(self):
        """returns a list of arrays
        each array contains all dots of with particular colour"""
        att = self._attributes
        att[np.where(att == None)] = "None"  # TODO check "is none"

        rtn = []
        for c in np.unique(att):
            if c is not None:
                da = DotArray(target_area_radius=self.target_area_radius,
                              min_dist_between=self.min_dist_between,
                              min_dist_area_boarder=self.min_dist_area_boarder)
                da.add(self.find(attribute=c))
                rtn.append(da)
        return rtn
