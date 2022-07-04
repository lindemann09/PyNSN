__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import copy
from random import shuffle

from ..exceptions import NoSolutionError, NoSizeDistributionError
from .. import _lib
from .size_distribution import SizeDistribution


#todo as_dict 
class NSNFactory(object):

    def __init__(self, target_area_radius, size_distribution=None,
                 min_dist_between=None,
                 min_dist_area_boarder=None):
        """

        Parameters
        ----------
        target_area_radius
        size_distribution
        min_dist_between
        min_dist_area_boarder
        """
        self._size_distr = None
        self.size_distribution = size_distribution # setter
        self._para = _lib.ArrayParameter(
            target_area_radius=target_area_radius,
            min_dist_between=min_dist_between,
            min_dist_area_boarder=min_dist_area_boarder)

    @property
    def array_parameter(self):
        return self._para

    @property
    def size_distribution(self):
        return self._size_distr

    @size_distribution.setter
    def size_distribution(self, sd_object):
        if sd_object is not None and not isinstance(sd_object,SizeDistribution):
            raise TypeError("size_distribution has to be of type SizeDistribution "
                            "or None, but not {}".format(
                type(sd_object).__name__))
        self._size_distr = sd_object

    def set_size_distribution(self, diameter=None, width=None, height=None,
                              rectangle_proportion=None):
        """

        Parameters
        ----------
        diameter
        width
        height
        rectangle_proportion

        Notes
        -----
        You can also also set size distribution by the property setter, which
        requires however as SizeDistribution object.
        """
        self._size_distr = SizeDistribution(
            diameter=diameter, width=width, height=height,
            rectangle_proportion=rectangle_proportion)

    def create_random_array(self, n_objects,
                            attributes=None,
                            allow_overlapping=False,
                            occupied_space=None):
        """
        occupied_space is a dot array (used for multicolour dot array (join after)

        attribute is an array, arrays are assigned randomly.


        Parameters
        ----------
        n_objects
        attributes
        allow_overlapping
        occupied_space

        Returns
        -------
        rtn : object array
        """
        if self._size_distr is None:
            raise NoSizeDistributionError("No size distribution defined. "
                                          "Use `set_size_distribution`.")
        if self._size_distr.diameter is not None:
            # DotArray
            rtn = _lib.DotArray(target_area_radius=self._para.target_area_radius,
                                min_dist_between=self._para.min_dist_between,
                                min_dist_area_boarder=self._para.min_dist_area_boarder)

            for dot in self._size_distr.sample(n=n_objects):
                try:
                    dot = rtn.get_random_free_position(ref_object=dot,
                                                       occupied_space=occupied_space,
                                                       allow_overlapping=allow_overlapping)
                except NoSolutionError as e:
                    raise NoSolutionError("Can't find a solution for {} items in this array".format(n_objects))
                rtn.add([dot])

        else:
            # RectArray
            rtn = _lib.RectangleArray(target_area_radius=self._para.target_area_radius,
                                      min_dist_between=self._para.min_dist_between,
                                      min_dist_area_boarder=self._para.min_dist_area_boarder)

            for rect in self._size_distr.sample(n=n_objects):
                try:
                    rect = rtn.get_random_free_position(ref_object=rect,
                                                        occupied_space=occupied_space,
                                                        allow_overlapping=allow_overlapping)
                except NoSolutionError:
                    raise NoSolutionError("Can't find a solution for {} ".format(n_objects) +
                                          "items in this array.")

                rtn.add([rect])

        # attribute assignment
        if isinstance(attributes, (tuple, list)):
            att = []
            while len(att) < n_objects:
                tmp = copy.copy(attributes)
                shuffle(tmp)
                att.extend(tmp)
            shuffle(att)
            rtn.set_attributes(att[:n_objects])
        else:
            rtn.set_attributes(attributes)

        return rtn

    def create_incremental_random_array(self, n_objects,
                                        attributes=None,
                                        allow_overlapping=False):
        """

        Parameters
        ----------
        n_objects
        attributes
        allow_overlapping

        Returns
        -------
        rtn : iterator of object arrays
        """
        previous = None
        for n in range(n_objects):
            current = self.create_random_array(n_objects=1,
                   attributes=attributes, allow_overlapping=allow_overlapping,
                   occupied_space=previous)
            if previous is not None:
                current.join(previous)
            previous = current
            yield current
