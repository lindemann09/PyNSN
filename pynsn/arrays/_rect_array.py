"""
Rectangle Array
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np

from ._generic_object_array import GenericObjectArray
from .._lib import misc
from ..shapes import Rectangle, Point
from . import _tools

class RectangleArray(GenericObjectArray):
    """
    """

    def __init__(self,
                 target_area_radius,
                 min_dist_between,
                 min_dist_area_boarder,
                 xy=None,
                 sizes=None,
                 attributes=None):
        """Rectangular array is restricted to a certain area, it has a target area
        and a minimum gap.

        This properties allows shuffling free position and adapting
        properties.

        """
        super().__init__(xy=xy, attributes=attributes,
                         target_area_radius=target_area_radius,
                         min_dist_between=min_dist_between,
                         min_dist_area_boarder=min_dist_area_boarder)
        self._sizes = np.array([])
        if sizes is not None:
            self._append_sizes(sizes)

        if self._xy.shape[0] != self._sizes.shape[0]:
            raise ValueError("Bad shaped data: " +
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
                  "min_dist_between": self.min_dist_between,
                  "target_area_radius": self.target_area_radius})
        return d

    def read_from_dict(self, the_dict):
        """read rectangle array from dict"""
        super().read_from_dict(the_dict)
        self._sizes = np.array(the_dict["sizes"])
        self.min_dist_between = the_dict["min_dist_between"]
        self.target_area_radius = the_dict["target_area_radius"]

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

        if self._properties.numerosity == 0:
            return RectangleArray(
                target_area_radius=self.target_area_radius,
                min_dist_between=self.min_dist_between)
        else:
            if indices is None:
                indices = list(range(self._properties.numerosity))

            return RectangleArray(
                        target_area_radius=self.target_area_radius,
                        min_dist_between=self.min_dist_between,
                        min_dist_area_boarder = self.min_dist_area_boarder,
                        xy=self._xy[indices, :].copy(),
                        sizes=self._sizes[indices].copy(),
                        attributes=self._attributes[indices].copy())

    def _xy_distances(self, rect):
        """return distances on both axes between rectangles and reference rec.
         negative number indicates overlap edges along that dimension.
        """
        if len(self._xy) == 0:
            return np.array([])
        else:
            pos_dist = np.abs(self._xy - rect.xy)
            max_overlap_dist = (self.sizes + rect.size)/2
            dist = pos_dist - max_overlap_dist
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
            # distances with a two one negative xy_distance --> overlap
            overlap = np.array(np.sum(d_xy<0, axis=1) == 2, dtype=int)
            return np.hypot(d_xy[:,0], d_xy[:,1])   * (overlap * -2 + 1) # overlap [0, 1, 0] -> [1,-1, 1]

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

    def find(self, size=None, attribute=None, edge=None):
        """returns indices of found objects

        2D-tuple
        """
        rtn = []
        for i in range(len(self._sizes)):
            if (size is not None and self._sizes[i] != size) or \
                    (attribute is not None and self._attributes[i] != attribute):
                continue
            rtn.append(i)

        if edge is not None:
            return rtn
        elif isinstance(edge, Point):
            new_rtn = []
            for i, rect in zip(rtn, self.get(indices=rtn)):
                if edge in list(rect.edges()):
                    new_rtn.append(i)
            return new_rtn
        else:
            raise TypeError("edge has to be of type Coordinate2D")

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
                rtn += "{},".format(self._properties.numerosity)
            rtn += "{},{},{},{}".format(self._xy[cnt, 0],
                                     self._xy[cnt, 1],
                                     self._sizes[cnt, 0],
                                     self._sizes[cnt, 1])
            if attribute_column:
                rtn += ", {}".format(self._attributes[cnt])
            rtn += "\n"
        return rtn

    def join(self, rect_array):
        """add another rect arrays"""
        assert isinstance(rect_array, RectangleArray)
        self.add(rect_array.get())

    def get_random_free_position(self,
                                 rectangle_size,
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
            convex_hull_xy = self.properties.convex_hull_positions.xy
        else:
            convex_hull_xy = None
        return _tools.get_random_free_position(
            the_object=Rectangle(xy=(0, 0), size=rectangle_size),
            target_area_radius = self.target_area_radius,
            allow_overlapping=allow_overlapping,
            distances_function=self.distances,
            min_dist_between=self.min_dist_between,
            min_dist_area_boarder=self.min_dist_area_boarder,
            occupied_space_distances_function=os_distance_fnc,
            convex_hull_xy=convex_hull_xy)

    def _remove_overlap_from_inner_to_outer(self):
        raise NotImplementedError()

    def realign(self):
        raise NotImplementedError()

    def shuffle_all_positions(self, allow_overlapping=False):
        raise NotImplementedError()

    def get_split_arrays(self):
        """returns a list of arrays
        each array contains all dots of with particular colour"""
        att = self._attributes
        att[np.where(att == None)] = "None"  # TODO check "is none"

        rtn = []
        for c in np.unique(att):
            if c is not None:
                da = RectangleArray(target_area_radius=self.target_area_radius,
                              min_dist_between=self.min_dist_between,
                              min_dist_area_boarder=self.min_dist_area_boarder)
                da.add(self.find(attribute=c))
                rtn.append(da)
        return rtn
