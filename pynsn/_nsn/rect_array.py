"""
Rectangle Array
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import random
import numpy as np
from scipy import spatial

from .cloud import _RestrictedCloud
from .._lib import misc
from .shape import Rectangle

class RectangleArray(_RestrictedCloud):
    """
    """

    def __init__(self, target_array_radius, minimum_gap,
                 xy=None, sizes=None, attributes=None):
        """Rectangular array is restricted to a certain area, it has a target area
        and a minimum gap.

        This features allows shuffling free position and matching
        features.

        """
        super().__init__(xy=xy, attributes=attributes,
                         target_array_radius=target_array_radius,
                         minimum_gap=minimum_gap)
        self._sizes = np.array([])
        if sizes is not None:
            self._append_sizes(sizes)

        if self._xy.shape[0] != self._sizes.shape[0]:
            raise RuntimeError("Bad shaped data: " +
                    u"xy has not the same length as sizes array")

    def _append_sizes(self, sizes):
        """returns number of added rows"""
        sizes = misc.numpy_array_2d(sizes)
        if len(self._sizes) == 0:
            empty = np.array([]).reshape((0, 2))  # ensure good shape of self.xy
            self._sizes = np.append(empty, sizes, axis=0)
        else:
            self._sizes = np.append(self._sizes, sizes, axis=0)
        return sizes.shape[0]

    def add(self, rectangles):
        """append one dot or list of dots"""
        if not isinstance(rectangles, (list, tuple)):
            rectangles = [rectangles]
        for r in rectangles:
            assert isinstance(r, Rectangle)
            self._append_xy_attribute(xy=r.xy, attributes=r.attribute)
            self._append_sizes((r.width, r.height))

    @property
    def sizes(self):
        return self._sizes

    @property
    def surface_areas(self):
        # a = w*h
        return self._sizes[:,0] * self._sizes[:,1]

    @property
    def perimeter(self):
        return 2 * (self._sizes[:,0] + self._sizes[:,1])

    def round(self, decimals=0, int_type=np.int32):
        """Round values of the array."""

        if decimals is None:
            return
        self._xy = misc.numpy_round2(self._xy, decimals=decimals,
                                     int_type=int_type)
        self._sizes = misc.numpy_round2(self._sizes, decimals=decimals,
                                     int_type=int_type)

    def as_dict(self):
        """
        """
        d = super().as_dict()
        d.update({"sizes": self._sizes.tolist(),
                  "minimum_gap": self.minimum_gap,
                  "target_array_radius": self.target_array_radius})
        return d

    def read_from_dict(self, the_dict):
        """read rectangle array from dict"""
        super().read_from_dict(the_dict)
        self._sizes = np.array(the_dict["sizes"])
        self.minimum_gap = the_dict["minimum_gap"]
        self.target_array_radius = the_dict["target_array_radius"]

    def clear(self):
        super().clear()
        self._sizes = np.array([])

    def delete(self, index):
        super().delete(index)
        self._sizes = np.delete(self._sizes, index)

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if indices is None:
            indices = list(range(self._features.numerosity))

        return RectangleArray(
                        target_array_radius=self.target_array_radius,
                        minimum_gap=self.minimum_gap,
                        xy=self._xy[indices, :].copy(),
                        sizes=self._sizes[indices].copy(),
                        attributes=self._attributes[indices].copy())

    def _xy_distances(self, rect):
        """return distances on both axes between rectangles and reference rec.
         0 indicates overlap off edges along that dimension.
        """
        if len(self._xy) == 0:
            return np.array([])
        else:
            pos_dist = self._xy - rect.xy
            max_overlap_dist = (self.sizes + rect.size)/2
            dist = pos_dist - max_overlap_dist
            dist[dist <  0] = 0
            return dist

    def distances(self, rect):
        """Euclidean Distances toward a single Rectangle
        negative numbers indicate overlap

        Returns
        -------
        distances : numpy array of distances
        """
        assert isinstance(rect, Rectangle)
        if len(self._xy) == 0:
            return np.array([])
        else:
            d_xy = self._xy_distances(rect)
            return np.hypot(d_xy[:,0], d_xy[:,1])

    @property
    def center_of_mass(self): # TODO does perimeter work for weighting?
        weighted_sum = np.sum(self._xy * self.perimeter[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self.perimeter)

    def get(self, indices=None):
        """returns all rectangles

         indices int or list of ints
         """

        if indices is None:
            return [Rectangle(xy=xy, size=s, attribute=att) \
                    for xy, s, att in zip(self._xy,
                                          self._sizes,
                                          self._attributes)]
        try:
            indices = list(indices)  # check if iterable
        except:
            indices = [indices]

        return [Rectangle(xy=xy, size=s, attribute=att) \
                    for xy, s, att in zip(self._xy[indices, :],
                                        self._sizes[indices],
                                        self._attributes[indices])]

    def find(self, size=None, attribute=None):
        """returns filtered dots
        """
        rtn = []
        for xy, s, att in zip(self._xy, self._sizes, self._attributes):
            if (size is not None and s != size) or \
                    (attribute is not None and att != attribute):
                continue

            rtn.append(Rectangle(xy=xy, size=s, attribute=att))
        return rtn

    def csv(self, variable_names=True,
            hash_column=True, num_idx_column=True,
            attribute_column=False):
        """Return the rectangle array as csv text

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
            rtn += u"x,y,width,height"
            if attribute_column:
                rtn += u",attribute"
            rtn += u"\n"

        obj_id = self.hash
        for cnt in range(len(self._xy)):
            if hash_column:
                rtn += "{0}, ".format(obj_id)
            if num_idx_column:
                rtn += "{},".format(self._features.numerosity)
            rtn += "{},{},{},{}".format(self._xy[cnt, 0],
                                     self._xy[cnt, 1],
                                     self._sizes[cnt, 0],
                                     self._sizes[cnt, 1])
            if attribute_column:
                rtn += ", {}".format(self._attributes[cnt])
            rtn += "\n"
        return rtn

    def join(self, rect_array):
        """add another dot arrays"""
        assert isinstance(rect_array, RectangleArray)
        self.add(rect_array.get())

    def random_free_position(self, rectangle_size,
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

        target_radius = self.target_array_radius - min_distance_area_boarder - \
                        min(rectangle_size)
        proposal_rect = Rectangle(xy=(0,0), size=rectangle_size)
        while True:
            cnt += 1
            ##  polar method seems to produce central clustering
            #  proposal_polar =  np.array([random.random(), random.random()]) *
            #                      (target_radius, TWO_PI)
            #proposal_xy = misc.polar2cartesian([proposal_polar])[0]
            #Note! np.random generates identical numbers under multiprocessing

            proposal_rect.xy = np.array([random.random(), random.random()]) \
                          * 2 * target_radius - target_radius

            bad_position = False
            if not squared_array:
                # check if one edge is outside
                for e in proposal_rect.edges():
                    if e.polar_radius >= target_radius:
                        bad_position = True
                        break

            if not bad_position and prefer_inside_field_area and \
                    cnt < try_out_inside_convex_hull:
                bad_position = delaunay.find_simplex(np.array(proposal_rect.xy)) < 0 # TODO check correctness, does it take into account size?

            if not bad_position and not allow_overlapping:
                # find bad_positions
                dist = self.distances(proposal_rect)
                if occupied_space:
                    dist = np.append(dist, occupied_space.distances(proposal_rect))
                idx = np.where(dist < self.minimum_gap)[0]  # overlapping dot ids
                bad_position = len(idx) > 0

            if not bad_position:
                return proposal_rect.xy
            elif cnt > 3000:
                raise RuntimeError(u"Can't find a free position") # TODO
