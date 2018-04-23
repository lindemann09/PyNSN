from __future__ import print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from . import random_beta
from copy import copy
from multiprocessing import Process, Event, Queue
from .dot_array import DotArray
from .item_features import ItemFeatures
from .dot_array_sequence import DASequence
from .log_file import GeneratorLogger
from . import continuous_property as cp


class DotArrayGenerator(object):

    def __init__(self,
                 max_array_radius,
                 dot_diameter_mean,
                 dot_diameter_range=None,
                 dot_diameter_std=None,
                 dot_colour=None,  # todo feature
                 minimum_gap=1,
                 **kwargs):

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

        if dot_diameter_std <= 0:
            dot_diameter_std = None
        if dot_diameter_range is not None and \
                (dot_diameter_mean <= dot_diameter_range[0] or
                         dot_diameter_mean >= dot_diameter_range[1] or
                         dot_diameter_range[0] >= dot_diameter_range[1]):
            raise RuntimeError("dot_diameter_mean has to be inside the defined dot_diameter_range")

        self.minimum_gap = minimum_gap
        self.max_array_radius = max_array_radius
        self.dot_diameter_range = dot_diameter_range
        self.dot_diameter_mean = dot_diameter_mean
        self.dot_diameter_std = dot_diameter_std
        self.item_feature = ItemFeatures(colour=dot_colour)

    def make(self, n_dots, occupied_space=None, logger=None):
        """occupied_space is a dot array (used for multicolour dot array (join after)"""

        rtn = DotArray(max_array_radius=self.max_array_radius,  # - distance_field_edge ?
                       minimum_gap=self.minimum_gap)

        for _ in range(n_dots):

            # diameter
            if self.dot_diameter_range is None or self.dot_diameter_std is None:
                # constant diameter
                diameter = self.dot_diameter_mean
            else:
                # draw diameter from beta distribution
                parameter = random_beta.shape_parameter_beta(self.dot_diameter_range,
                                                             self.dot_diameter_mean,
                                                             self.dot_diameter_std)
                diameter = random_beta.random_beta(
                    self.dot_diameter_range, parameter)

            xy = rtn.random_free_dot_position(dot_diameter=diameter, occupied_space=occupied_space)
            rtn.append(xy=xy, diameters=diameter, features=self.item_feature)

        if logger is not None:
            if not isinstance(logger, GeneratorLogger):
                raise RuntimeError("logger has to be None or a GeneratorLogger")
            logger.log(rtn)

        return rtn

    def as_dict(self):
        return {"max_array_radius": self.max_array_radius,
             "dot_diameter_mean": self.dot_diameter_mean,
             "dot_diameter_range": self.dot_diameter_range,
             "dot_diameter_std": self.dot_diameter_std,
             "dot_colour": self.item_feature.colour.colour,  ##todo feature
             "minimum_gap": self.minimum_gap}


def make_dot_array_sequence(reference_dot_array,
                            match_properties,
                            min_max_numerosity,
                            extra_space,  # fitting convex hull and density might result in enlarged arrays
                            center_array=True,
                            logger=None):
    """Methods takes take , you might use make Process
        match_properties:
                continuous property or list of continuous properties to be match
                or None
     returns False is error occured (see self.error)
    """

    if not isinstance(reference_dot_array, DotArray):
        raise TypeError("Reference_dot_array has to be DotArray, but not {}".format(
            type(reference_dot_array).__name__))

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

    cp.check_list_continuous_properties(match_properties, check_set_value=False)

    # copy and change values to match this stimulus
    match_props = []
    for m in match_properties:
        m = copy(m)
        m.set_value(reference_dot_array)
        match_props.append(m)

    # adjust reference (basically centering)
    reference_da = reference_dot_array.copy()
    reference_da.max_array_radius += (extra_space // 2)  # add extra space
    if center_array:
        reference_da._xy -= reference_da.center_of_outer_positions
        reference_da.set_array_modified()

    # matched deviants
    prefer_keeping_convex_hull = False
    for x in match_properties:
        if isinstance(x, cp.ConvexHull) or \
                (isinstance(x, cp.Density) and x.match_ratio_convhull2area < 1):
            prefer_keeping_convex_hull = True
            break

    rtn = DASequence()
    rtn.method = match_props

    min, max = sorted(min_max_numerosity)
    # decreasing
    if min < reference_dot_array.prop_numerosity:
        da_sequence, error = _make_matched_deviants(reference_da=reference_da,
                                                    match_props=match_props,
                                                    target_numerosity=min,
                                                    prefer_keeping_convex_hull=prefer_keeping_convex_hull)
        rtn.append_dot_arrays(list(reversed(da_sequence)))
        if error is not None:
            rtn.error = error
    # reference
    rtn.append_dot_arrays(reference_da)
    # increasing
    if max > reference_dot_array.prop_numerosity:
        da_sequence, error = _make_matched_deviants(reference_da=reference_da,
                                                    match_props=match_props,
                                                    target_numerosity=max,
                                                    prefer_keeping_convex_hull=prefer_keeping_convex_hull)
        rtn.append_dot_arrays(da_sequence)
        if error is not None:
            rtn.error = error

    if logger is not None:
        if not isinstance(logger, GeneratorLogger):
            raise TypeError("logger has to be None or a GeneratorLogger, and not {}".format(
                type(logger).__name__))

        logger.log(rtn)

    return rtn


def _make_matched_deviants(reference_da,
                           match_props,
                           target_numerosity,
                           prefer_keeping_convex_hull):  # TODO center array OK?
    """helper function. Do not use this method. Please use make"""

    if reference_da.prop_numerosity == target_numerosity:
        change = 0
    elif reference_da.prop_numerosity > target_numerosity:
        change = -1
    else:
        change = 1

    da = reference_da.copy()
    da_sequence = []

    error = None
    while True:
        da = da.number_deviant(change_numerosity=change,
                               prefer_keeping_convex_hull=prefer_keeping_convex_hull)

        if len(match_props) > 0:
            da.match(match_props, center_array=False)

        cnt = 0
        while True:
            cnt += 1
            ok, mesg = da.realign()
            if ok:
                break
            if cnt > 10:
                error = u"ERROR: realign, " + str(cnt) + ", " + str(da.prop_numerosity)

        da_sequence.append(da)

        if error is not None or da.prop_numerosity == target_numerosity:
            break

    return da_sequence, error


class DASequenceMakeProcess(Process):
    def __init__(self, reference_dot_array, min_max_numerosity, match_properties, extra_space,
                 n_trials=3, logger=None):

        """
        property: da_sequence, after processes finished
        Event(): data_available
        """

        super(DASequenceMakeProcess, self).__init__()

        self.data_available = Event()
        self._data_queue = Queue()
        self._da_sequence = None

        if not isinstance(match_properties, (tuple, list)):
            match_properties = [match_properties]
        cp.check_list_continuous_properties(match_properties)
        self.match_properties = match_properties
        self.reference_dot_array = reference_dot_array
        self.min_max_numerosity = min_max_numerosity
        self.extra_space = extra_space

        if n_trials < 1:
            self._n_tryouts = 1
        else:
            self._n_tryouts = n_trials

        self.logger = logger
        if not isinstance(logger, (type(None), GeneratorLogger)):
            raise RuntimeError(u"logger has to be None or a GeneratorLogger")

    @property
    def da_sequence(self):
        self._read_queue()
        return self._da_sequence

    def join(self, timeout=1):
        self._read_queue()
        super(DASequenceMakeProcess, self).join(timeout)

    def _read_queue(self):
        if self._da_sequence is not None:
            return
        if self.is_alive() or not self._data_queue.empty():
            self._da_sequence = self._data_queue.get()

    def run(self):
        cnt = 0
        da_seq = None
        while cnt < self._n_tryouts:
            cnt += 1
            da_seq = make_dot_array_sequence(reference_dot_array=self.reference_dot_array,
                                             match_properties=self.match_properties,
                                             extra_space=self.extra_space,
                                             min_max_numerosity=self.min_max_numerosity,
                                             logger=self.logger)
            if da_seq.error is None:
                break
                # print("remix")

        self.data_available.set()
        self._data_queue.put(da_seq)
