__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import deepcopy
import random
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import spatial

from . import geometry
from . import shapes
from .misc import NoSolutionError
from .base_array import BaseArray
from ..visual_properties import fit


class ABCObjectArray(BaseArray, metaclass=ABCMeta):

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
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """
        pass

    @abstractmethod
    def iter_objects(self, indices=None):
        pass

    @abstractmethod
    def add(self, something):
        pass

    @abstractmethod
    def find(self, attribute):
        pass

    @abstractmethod
    def distances(self, ref_object):
        # override this method
        pass

    @abstractmethod
    def check_stand_outs(self):
        pass

    def join(self, object_array):
        """add another object arrays"""
        self.add(object_array.iter_objects())

    def center_of_field_area(self):
        return geometry.center_of_positions(self.properties.convex_hull.xy)

    def distances_matrix(self, between_positions=False):
        """between position ignores the dot size"""
        if between_positions:
            return spatial.distance.cdist(self._xy, self._xy)
        # matrix with all distance between all points
        dist = np.array([self.distances(d) for d in self.iter_objects()])
        return dist

    def check_overlaps(self):
        """return pairs of indices of overlapping of objects
        numpy.array
        """
        dist = self.distances_matrix(between_positions=False)
        overlap = np.where(np.triu(dist, k=1) < 0)
        return np.array(overlap).T

    def center_of_mass(self):
        weighted_sum = np.sum(self._xy * self.perimeter[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self.perimeter)

    def get_random_free_position(self, ref_object,
                                 allow_overlapping=False,
                                 inside_convex_hull=False,
                                 occupied_space=None):
        """returns the copy of object of at a random free position

        raise exception if not found
        occupied space: see generator generate
        """

        N_ATTEMPTS = 3000

        if isinstance(ref_object, shapes.Dot):
            object_size = ref_object.diameter / 2.0
        elif isinstance(ref_object, shapes.Rectangle):
            object_size = max(ref_object.size)
        else:
            raise NotImplementedError("Not implemented for {}".format(
                type(ref_object).__name__))
        if occupied_space is not None and \
                not isinstance(occupied_space, ABCObjectArray): #FIXME check
            raise TypeError("Occupied_space has to be a Dot or Rectangle Array or None.")

        area_rad = self.target_area_radius - self.min_dist_area_boarder - object_size
        rtn_object = deepcopy(ref_object)  # tested deepcopy required

        cnt = 0
        while True:
            cnt += 1
            ##  polar method seems to produce central clustering
            #  proposal_polar =  np.array([random.random(), random.random()]) *
            #                      (target_radius, TWO_PI)
            rtn_object.xy = np.array([random.random(), random.random()]) \
                            * 2 * area_rad - area_rad

            # is outside area
            if isinstance(ref_object, shapes.Dot):
                bad_position = area_rad <= rtn_object.polar_radius
            else:
                # Rect: check if one edge is outside
                bad_position = False
                for e in rtn_object.iter_edges():
                    if e.polar_radius >= area_rad:
                        bad_position = True
                        break

            if not bad_position and not allow_overlapping:
                # find bad_positions
                dist = self.distances(rtn_object)
                if isinstance(occupied_space, ABCObjectArray):
                    dist = np.append(dist, occupied_space.distances(rtn_object))
                bad_position = sum(dist < self.min_dist_between) > 0  # at least one is overlapping

            if not bad_position and inside_convex_hull:
                # use only those that do not change the convex hull
                tmp_array = self.copy(deepcopy=True)
                tmp_array.add([rtn_object])
                bad_position = tmp_array.properties.convex_hull != \
                               self.properties.convex_hull

            if not bad_position:
                return rtn_object
            elif cnt > N_ATTEMPTS:
                raise NoSolutionError(u"Can't find a free position")

    def shuffle_all_positions(self, allow_overlapping=False):
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        all_objects = list(self.iter_objects())
        self.clear()
        for obj in all_objects:
            try:
                new = self.get_random_free_position(obj,
                                        allow_overlapping=allow_overlapping)
            except NoSolutionError as e:
                raise NoSolutionError("Can't shuffle dot array. No free positions found.")
            self.add([new])


    def get_number_deviant(self, change_numerosity, preserve_field_area=False):
        """number deviant
        """
        object_array = self.copy()
        new_num = self.properties.numerosity + change_numerosity
        fit.numerosity(object_array, value=new_num,
                       center_of_field_area=preserve_field_area)
        return object_array

    def replace_overlapping_objects(self, center_of_field_area=False,
                                    lenient=True):
        """
        Returns
        Parameters
        ----------
        center_of_field_area
        lenient

        Returns
        -------
        convex_hull_had_changed
        """

        warning_info = "Field area had to change, because two overlapping " + \
                       "objects on convex hull"
        convex_hull_had_changed = False

        overlaps = self.check_overlaps()
        while len(overlaps):
            if center_of_field_area:
                # check if overlaps are on convex hull
                # do not replace convexhull objects and try to
                # take that one not on convex hull or raise error/warning
                ch_idx = self.properties.convex_hull.indices
                if overlaps[0, 0] not in ch_idx:
                    idx = overlaps[0, 0]
                elif overlaps[0, 1] not in ch_idx:
                    idx = overlaps[0, 1]
                elif lenient:
                    # warning
                    convex_hull_had_changed = True
                    idx = overlaps[0, 0]
                else:
                    raise NoSolutionError("Can't replace overlap and keep convex hull unchanged. " +
                                          warning_info)
            else:
                idx = overlaps[0, 0]

            obj = next(self.iter_objects(idx))
            self.delete(idx)
            obj = self.get_random_free_position(ref_object=obj,
                                                inside_convex_hull=center_of_field_area)
            self.add([obj])
            overlaps = self.check_overlaps()

        if convex_hull_had_changed:
            print("Warning: " + warning_info)

        return convex_hull_had_changed

    @abstractmethod
    def center_array(self):
        """places array in target area as central and possible and tries to
        remove any stand_outs"""
        pass

    @abstractmethod
    def realign(self):
        raise NotImplementedError()

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
                da.add(self.find(attribute=c))
                rtn.append(da)
        return rtn


# TODO  everywhere: file header doc and author information
