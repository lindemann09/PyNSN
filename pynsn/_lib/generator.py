from __future__ import print_function, division
from builtins import map, filter, range, zip

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from . import misc
from copy import copy
from multiprocessing import Pool
from .dot_array import DotArray
from .item_attributes import ItemAttributes
from .dot_array_sequence import DASequence
from .log_file import GeneratorLogger
from . import features as cp


class DotArrayGenerator(object):

    def __init__(self,
                 target_area_radius,
                 item_diameter_mean,
                 item_diameter_range=None,
                 item_diameter_std=None,
                 item_colour=None,  # todo feature
                 minimum_gap=1): # TODO check minim gap

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

    def make(self, n_dots, occupied_space=None, logger=None):
        """occupied_space is a dot array (used for multicolour dot array (join after)"""

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

            xy = rtn.random_free_dot_position(dot_diameter=diameter, occupied_space=occupied_space)
            rtn.append(xy=xy, item_diameters=diameter, features=self.item_feature)

        if logger is not None:
            if not isinstance(logger, GeneratorLogger):
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


    def make_iter(self, list_of_n_dots, occupied_space=None, logger=None, multiprocessing=False):  # TODO  never checked
        args = map(lambda x: (self, x, occupied_space, logger), list_of_n_dots)

        if multiprocessing:
            return Pool().imap(_make_imap_helper, args)
        else:
            return map(_make_imap_helper, args)


def _make_imap_helper(args):
    generator = args[0]
    return generator.make(reference_dot_array=args[1],
                          occupied_space=args[2],
                          logger=args[3])


class DASequenceGenerator(object):

    def __init__(self,
                 match_properties,
                 min_max_numerosity,
                 extra_space,  # fitting convex hull and density might result in enlarged arrays
                 center_array=True):
        """Methods takes take , you might use make Process
            match_properties:
                    continuous property or list of continuous properties to be match
                    or None
         returns False is error occured (see self.error)
        """

        try:
            l = len(min_max_numerosity)
        except:
            l = 0
        if l != 2:
            raise ValueError("min_max_numerosity has to be a pair of (min, max)")

        if match_properties is None:
            match_properties = []
        elif not isinstance(match_properties, (tuple, list)):
            match_properties = [match_properties]
        cp.check_feature_list(match_properties, check_set_value=False)

        self.match_properties = match_properties
        self.min_max_numerosity = min_max_numerosity
        self.extra_space = extra_space
        self.center_array = center_array


    def make(self, reference_dot_array, logger=None): # todo could be an iterator
        """Methods takes take , you might use make Process
            match_properties:
                    continuous property or list of continuous properties to be match
                    or None
         returns False is error occured (see self.error)
        """

        if not isinstance(reference_dot_array, DotArray):
            raise TypeError("Reference_dot_array has to be DotArray, but not {}".format(
                type(reference_dot_array).__name__))

        # copy and change values to match this stimulus
        match_props = []
        prefer_keeping_field_area = False
        for m in self.match_properties:
            m = copy(m)
            m.set_value(reference_dot_array)
            match_props.append(m)
            if isinstance(m, cp.LogSpacing().dependencies ) or \
                    (isinstance(m, cp.Coverage) and m.match_ratio_fieldarea2totalarea < 1):
                prefer_keeping_field_area = True
                break

        # adjust reference (basically centering)
        reference_da = reference_dot_array.copy()
        reference_da.target_array_radius += (self.extra_space // 2)  # add extra space
        if self.center_array:
            reference_da._xy -= reference_da.center_of_outer_positions
            reference_da.set_array_modified()

        # matched deviants
        rtn = DASequence()
        rtn.method = match_props

        min, max = sorted(self.min_max_numerosity)
        # decreasing
        if min < reference_dot_array.feature_numerosity:
            da_sequence, error = DASequenceGenerator._make_matched_deviants(
                reference_da=reference_da,
                match_props=match_props,
                target_numerosity=min,
                prefer_keeping_field_area=prefer_keeping_field_area)
            rtn.append_dot_arrays(list(reversed(da_sequence)))
            if error is not None:
                rtn.error = error
        # reference
        rtn.append_dot_arrays(reference_da)
        # increasing
        if max > reference_dot_array.feature_numerosity:
            da_sequence, error = DASequenceGenerator._make_matched_deviants(
                reference_da=reference_da,
                match_props=match_props,
                target_numerosity=max,
                prefer_keeping_field_area=prefer_keeping_field_area)
            rtn.append_dot_arrays(da_sequence)
            if error is not None:
                rtn.error = error

        if logger is not None:
            if not isinstance(logger, GeneratorLogger):
                raise TypeError("logger has to be None or a GeneratorLogger, and not {}".format(
                    type(logger).__name__))

            logger.log(rtn)

        return rtn

    @staticmethod
    def _make_matched_deviants(reference_da, match_props, target_numerosity,
                               prefer_keeping_field_area):  # TODO center array OK?
        """helper function. Do not use this method. Please use make"""

        if reference_da.feature_numerosity == target_numerosity:
            change = 0
        elif reference_da.feature_numerosity > target_numerosity:
            change = -1
        else:
            change = 1

        da = reference_da.copy()
        da_sequence = []

        error = None
        while True:
            da = da.number_deviant(change_numerosity=change,
                                   prefer_keeping_field_area=prefer_keeping_field_area)

            if len(match_props) > 0:
                da.match(match_props, center_array=False)

            cnt = 0
            while True:
                cnt += 1
                ok, mesg = da.realign()
                if ok:
                    break
                if cnt > 10:
                    error = u"ERROR: realign, " + str(cnt) + ", " + str(da.feature_numerosity)

            da_sequence.append(da)

            if error is not None or da.feature_numerosity == target_numerosity:
                break

        return da_sequence, error


# class DASequenceMakeProcess(Process):
#     def __init__(self, reference_dot_array, min_max_numerosity, match_properties, extra_space,
#                  n_trials=3, logger=None):
#
#         """
#         property: da_sequence, after processes finished
#         Event(): data_available
#         """
#
#         super(DASequenceMakeProcess, self).__init__()
#
#         self.data_available = Event()
#         self._data_queue = Queue()
#         self._da_sequence = None
#
#         if not isinstance(match_properties, (tuple, list)):
#             match_properties = [match_properties]
#         cp.check_list_continuous_properties(match_properties)
#         self.match_properties = match_properties
#         self.reference_dot_array = reference_dot_array
#         self.min_max_numerosity = min_max_numerosity
#         self.extra_space = extra_space
#
#         if n_trials < 1:
#             self._n_tryouts = 1
#         else:
#             self._n_tryouts = n_trials
#
#         self.logger = logger
#         if not isinstance(logger, (type(None), GeneratorLogger)):
#             raise RuntimeError(u"logger has to be None or a GeneratorLogger")
#
#     @property
#     def da_sequence(self):
#         self._read_queue()
#         return self._da_sequence
#
#     def join(self, timeout=1):
#         self._read_queue()
#         super(DASequenceMakeProcess, self).join(timeout)
#
#     def _read_queue(self):
#         if self._da_sequence is not None:
#             return
#         if self.is_alive() or not self._data_queue.empty():
#             self._da_sequence = self._data_queue.get()
#
#     def run(self):
#         cnt = 0
#         da_seq = None
#         while cnt < self._n_tryouts:
#             cnt += 1
#             da_seq = make_dot_array_sequence(reference_dot_array=self.reference_dot_array,
#                                              match_properties=self.match_properties,
#                                              extra_space=self.extra_space,
#                                              min_max_numerosity=self.min_max_numerosity,
#                                              logger=self.logger)
#             if da_seq.error is None:
#                 break
#                 # print("remix")
#
#         self.data_available.set()
#         self._data_queue.put(da_seq)
#

# class GeneratorProcess(Process): # TODO required, better get building iterators is task_nsn_production
#     # calls the make method of the generator
#     def __init__(self, generator, **kwargs):
#         Process.__init__(self)
#         self._generator = generator
#         self._kwargs = kwargs
#         self._queue = Queue()
#         self._result = None
#         self.data_available = Event()
#         self.deamon = True
#
#     @property
#     def result(self):
#         if self._result is not None:
#             return self._result
#         else:
#             try:
#                 self._result = self._queue.get_nowait()
#             except:
#                 pass
#         return self._result
#
#     def run(self):
#         # print(self)
#         rtn = self._generator.make(**self._kwargs)
#         self.data_available.set()
#         self._queue.put_nowait(rtn)
#
#     def join(self, timeout=1):
#
#         # print("wait join :" + str(self))
#         self._result = self.result  # ensure ready result
#         while True:
#             try:  # empty queque
#                 self._queue.get_nowait()
#             except:
#                 break
#         self.terminate()
#         Process.join(self, timeout=timeout)
#         # print("joined:" + str(self))
