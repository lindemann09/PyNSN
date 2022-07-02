"""
Rectangle Array
"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np

from ._base_array import ABCObjectArray
from .._lib import misc
from ..shapes import Rectangle, Point


class RectangleArray(ABCObjectArray):
    """
    """

    def __init__(self,
                 target_area_radius,
                 min_dist_between=2,
                 min_dist_area_boarder=1,
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
        try:
            rectangles = list(rectangles)
        except TypeError:
            rectangles = [rectangles]

        for r in rectangles:
            assert isinstance(r, Rectangle)
            self._append_xy_attribute(xy=r.xy, attributes=r.attribute)
            self._append_sizes((r.width, r.height))
        self.properties.reset()

    @property
    def sizes(self):
        return self._sizes

    @property
    def surface_areas(self):
        # a = w*h
        return self._sizes[:, 0] * self._sizes[:, 1]

    @property
    def perimeter(self):
        return 2 * (self._sizes[:, 0] + self._sizes[:, 1])

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
        d.update({"sizes": self._sizes.tolist()})
        return d

    def read_from_dict(self, the_dict):
        """read rectangle array from dict"""
        super().read_from_dict(the_dict)
        self._sizes = np.array(the_dict["sizes"])

    def clear(self):
        super().clear()
        self._sizes = np.array([])

    def delete(self, index):
        super().delete(index)
        self._sizes = np.delete(self._sizes, index)

    def copy(self, indices=None, deepcopy=True):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._xy) == 0:
            return RectangleArray(
                target_area_radius=self.target_area_radius,
                min_dist_area_boarder=self.min_dist_area_boarder,
                min_dist_between=self.min_dist_between)
        if indices is None:
            indices = list(range(len(self._xy)))

        if deepcopy:
            return RectangleArray(
                target_area_radius=self.target_area_radius,
                min_dist_between=self.min_dist_between,
                min_dist_area_boarder=self.min_dist_area_boarder,
                xy=self._xy[indices, :].copy(),
                sizes=self._sizes[indices].copy(),
                attributes=self._attributes[indices].copy())
        else:
            return RectangleArray(
                target_area_radius=self.target_area_radius,
                min_dist_between=self.min_dist_between,
                min_dist_area_boarder=self.min_dist_area_boarder,
                xy=self._xy[indices, :],
                sizes=self._sizes[indices],
                attributes=self._attributes[indices])

    def _xy_distances(self, rect):
        """return distances on both axes between rectangles and reference rec.
         negative number indicates overlap edges along that dimension.
        """
        if len(self._xy) == 0:
            return np.array([])
        else:
            pos_dist = np.abs(self._xy - rect.xy)
            max_overlap_dist = (self.sizes + rect.size) / 2
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
            # distances with a two one negative xy_distances --> overlap
            overlap = np.array(np.sum(d_xy < 0, axis=1) == 2, dtype=int)
            return np.hypot(d_xy[:, 0], d_xy[:, 1]) * (overlap * -2 + 1)  # overlap [0, 1, 0] -> [1,-1, 1]

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
            data = zip(self._xy, self._sizes, self._attributes)
        else:
            try:
                indices = list(indices)  # check if iterable
            except TypeError:
                indices = [indices]
            data = zip(self._xy[indices, :], self._sizes[indices],
                       self._attributes[indices])

        for xy, s, att in data:
            rtn = Rectangle(xy=xy, size=s, attribute=att)
            yield rtn

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
            for i, rect in zip(rtn, self.iter_objects(indices=rtn)):
                if edge in list(rect.iter_edges()):
                    new_rtn.append(i)
            return new_rtn
        else:
            raise TypeError("edge has to be of type Points")

    def csv(self, variable_names=True,
            hash_column=True,
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
                rtn += "hash,"
            rtn += "x,y,width,height"
            if attribute_column:
                rtn += ",attribute"
            rtn += "\n"

        obj_id = self.hash
        for cnt in range(len(self._xy)):
            if hash_column:
                rtn += "{0}, ".format(obj_id)
            rtn += "{},{},{},{}".format(self._xy[cnt, 0],
                                        self._xy[cnt, 1],
                                        self._sizes[cnt, 0],
                                        self._sizes[cnt, 1])
            if attribute_column:
                rtn += ", {}".format(self._attributes[cnt])
            rtn += "\n"
        return rtn

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

    def center_array(self):
        raise NotImplementedError()

    def check_stand_outs(self):
        raise NotImplementedError()

    def realign(self):
        raise NotImplementedError()

# FIXME overlapping does not work for rectangles
# FIXME picture random size distribution