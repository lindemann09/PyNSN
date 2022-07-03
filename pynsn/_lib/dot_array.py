"""
Dot Array
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np

from .object_array import ABCObjectArray
from .._lib import misc, geometry
from .shapes import Dot
from . import array_tools

# TODO: How to deal with rounding? Is saving to precises? Suggestion:
#  introduction precision parameter that is used by as_dict and get_csv and
#  hash


class DotArray(ABCObjectArray):
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

        This properties allows shuffling free position and adapting
        properties.
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
        try:
            dots = list(dots)
        except TypeError:
            dots = [dots]
        for d in dots:
            assert isinstance(d, Dot)
            self._append_xy_attribute(xy=d.xy, attributes=d.attribute)
            self._diameters = np.append(self._diameters, d.diameter)
        self.properties.reset()

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
        d.update({"diameters": self._diameters.tolist()})
        return d

    def read_from_dict(self, the_dict):
        """read Dot collection from dict"""
        super().read_from_dict(the_dict)
        self._diameters = np.array(the_dict["diameters"])

    def clear(self):
        super().clear()
        self._diameters = np.array([])

    def delete(self, index):
        super().delete(index)
        self._diameters = np.delete(self._diameters, index)

    def copy(self, indices=None, deepcopy=True):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._xy) == 0:
            return DotArray(target_area_radius=self.target_area_radius,
                            min_dist_between=self.min_dist_between,
                            min_dist_area_boarder=self.min_dist_area_boarder)

        if indices is None:
            indices = list(range(len(self._xy) ))

        if deepcopy:
            return DotArray(target_area_radius=self.target_area_radius,
                        min_dist_between=self.min_dist_between,
                        min_dist_area_boarder = self.min_dist_area_boarder,
                        xy=self._xy[indices, :].copy(),
                        diameters=self._diameters[indices].copy(),
                        attributes=self._attributes[indices].copy())
        else:
            return DotArray(target_area_radius=self.target_area_radius,
                        min_dist_between=self.min_dist_between,
                        min_dist_area_boarder = self.min_dist_area_boarder,
                        xy=self._xy[indices, :],
                        diameters=self._diameters[indices],
                        attributes=self._attributes[indices])

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

    def iter_objects(self, indices=None):
        """iterate over all or a part of the objects

        Parameters
        ----------
        indices

        Notes
        -----
        To iterate all object you might all use the class iterator __iter__:
        >>> for obj in my_array:
        >>>    print(obj)
        """

        if indices is None:
            data = zip(self._xy, self._diameters, self._attributes)
        else:
            try:
                indices = list(indices)  # check if iterable
            except TypeError:
                indices = [indices]
            data = zip(self._xy[indices, :], self._diameters[indices],
                       self._attributes[indices])

        for xy, dia, att in data:
            yield Dot(xy=xy, diameter=dia, attribute=att)

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

    def csv(self, variable_names=True, hash_column=False,
            attribute_column=False):
        """Return the dot array as csv text

        Parameter
        ---------
        variable_names : bool, optional
            if True variable name will be printed in the first line
        """
        size_dict = {"diameter": self._diameters}
        if attribute_column:
            attr = self.attributes
        else:
            attr = None
        if hash_column:
            array_hash = self.hash
        else:
            array_hash = None
        return misc.make_csv(xy=self._xy,
                             size_data_dict=size_dict,
                             attributes=attr, array_hash=array_hash,
                             make_variable_names=variable_names)

    def __realign_old(self):
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

        xy, shift_required = array_tools.remove_overlap_from_inner_to_outer(
            xy=self.xy, min_dist_between=self.min_dist_between,
            distance_matrix_function=self.distances_matrix)

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
                   str(self._properties.numerosity) + ")"

        self._properties.reset()
        if not shift_required:
            return True, ""
        else:
            return self.realign()  # recursion

    def check_stand_outs(self):  # FIXME ignores object size
        """returns indices of object that stand out"""
        ch_radii = geometry.cartesian2polar(self._xy, radii_only=True)
        return np.where(ch_radii > self.target_area_radius -
                        self.min_dist_area_boarder)[0]

    def center_array(self):
        """places array in target area as central and possible and tries to
        remove any stand_outs"""
        raise NotImplementedError()

    def realign(self):
        error = False
        realign_required = False

        self.center_array()
        cnt = 0
        while True:
            cnt += 1
            if cnt > 20:
                error = True
                break

        if error:
            raise RuntimeError("Can't find solution when removing outlier (n=" + \
                   str(self._properties.numerosity) + ")")

        self._properties.reset()
        if not realign_required:
            return True, ""
        else:
            return self.realign()  # recursion