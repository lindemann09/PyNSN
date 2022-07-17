__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from ..exceptions import NoSolutionError, NoAppearanceDefinedError
from .base_classes import ArrayParameter
from .._lib.misc import dict_to_text
from .dot_array import DotArray
from .rect_array import RectangleArray
from .appearance_sampler import AppearanceSampler


class NSNFactory(ArrayParameter, AppearanceSampler):

    def __init__(self, target_area_radius,
                 min_dist_between=None,
                 min_dist_area_boarder=None):
        """

        Parameters
        ----------
        target_area_radius
        min_dist_between
        min_dist_area_boarder
        """

        AppearanceSampler.__init__(self)
        ArrayParameter.__init__(self,
            target_area_radius=target_area_radius,
            min_dist_between=min_dist_between,
            min_dist_area_boarder=min_dist_area_boarder)

    def create_random_array(self, n_objects,
                            allow_overlapping=False,
                            occupied_space=None):
        """
        occupied_space is a dot array (used for multicolour dot array (join after)

        attribute is an array, arrays are assigned randomly.


        Parameters
        ----------
        n_objects
        allow_overlapping
        occupied_space

        Returns
        -------
        rtn : object array
        """
        if not self.is_appearance_set:
            raise NoAppearanceDefinedError("No appearance defined. Please use 'set_dot', 'set_rect_or' "
                                           "'set_appearance'")
        if self._diameter is not None:
            # DotArray
            rtn = DotArray(target_area_radius=self.target_area_radius,
                                min_dist_between=self.min_dist_between,
                                min_dist_area_boarder=self.min_dist_area_boarder)

            for dot in self.sample(n=n_objects):
                try:
                    dot = rtn.get_random_free_position(ref_object=dot,
                                                       occupied_space=occupied_space,
                                                       allow_overlapping=allow_overlapping)
                except NoSolutionError as e:
                    raise NoSolutionError("Can't find a solution for {} items in this array".format(n_objects))
                rtn.add([dot])

        else:
            # RectArray
            rtn = RectangleArray(target_area_radius=self.target_area_radius,
                                      min_dist_between=self.min_dist_between,
                                      min_dist_area_boarder=self.min_dist_area_boarder)

            for rect in self.sample(n=n_objects):
                try:
                    rect = rtn.get_random_free_position(ref_object=rect,
                                                        occupied_space=occupied_space,
                                                        allow_overlapping=allow_overlapping)
                except NoSolutionError:
                    raise NoSolutionError("Can't find a solution for {} ".format(n_objects) +
                                          "items in this array.")

                rtn.add([rect])
        return rtn

    def create_incremental_random_array(self, n_objects,
                                        allow_overlapping=False):
        """

        Parameters
        ----------
        n_objects
        allow_overlapping

        Returns
        -------
        rtn : iterator of object arrays
        """
        previous = None
        for n in range(n_objects):
            current = self.create_random_array(n_objects=1,
                   allow_overlapping=allow_overlapping,
                   occupied_space=previous)
            if previous is not None:
                current.join(previous)
            previous = current
            yield current

    def as_dict(self):
        d = ArrayParameter.as_dict(self)
        d.update(AppearanceSampler.as_dict(self))
        return d

    def __str__(self):
        return dict_to_text(self.as_dict(), col_a=12, col_b=7)