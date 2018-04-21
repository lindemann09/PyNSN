from __future__ import  print_function, division
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
                 dot_colour = None,
                 minimum_gap=1,
                 logger=None):

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
        self.set_logger(logger)


    def set_logger(self, logger):
        self.logger = logger
        if not isinstance(logger, (type(None), GeneratorLogger)):
            raise RuntimeError("logger has to be None or a GeneratorLogger")


    def make(self, n_dots, inhibit_logging=False):

        rtn = DotArray(max_array_radius= self.max_array_radius,  # - distance_field_edge ?
                       minimum_gap=self.minimum_gap)

        for _ in range(n_dots):

            #diameter
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

            xy = rtn.random_free_dot_position(dot_diameter=diameter)
            rtn.append(xy=xy, diameters=diameter, features=self.item_feature)

        if not inhibit_logging and self.logger is not None:
            self.logger.log(rtn)

        return rtn




class DASequenceGenerator(object):

    def __init__(self, max_dot_array, logger=None):
        """  makes sequence of deviants by subtracting dots

            sqeeze factor: when adapting for convex hull, few point shift excentrically, it is
                  therefore usefull to sqeeze the stimulus before. We do that therefore all stimuli

        Stimulus will be center before making variants

        """
        self.max_dot_array=max_dot_array
        self.logger = logger

    @property
    def logger(self):
        return self._logger

    @logger.setter
    def logger(self, x):
        if not isinstance(x, (type(None), GeneratorLogger)):
            raise TypeError("logger has to be None or a GeneratorLogger, and not {}".format(
                type(x).__name__))
        self._logger = x

    @property
    def max_dot_array(self):
        return self._da

    @max_dot_array.setter
    def max_dot_array(self, x):
        if not isinstance(x, DotArray):
            raise TypeError("Max_dot_array has to be DotArray, but not {}".format(
                            type(x).__name__))
        self._da = x.copy()


    def make(self, match_properties, min_numerosity,
             extra_space,  # fitting convex hull and density might result in enlarged arrays
             inhibit_logging=False):
        """Methods takes take , you might use make Process
            match_properties:
                    continuous property or list of continuous properties to be match
                    or None
         returns False is error occured (see self.error)
        """

        if match_properties is None:
            match_properties = []
        elif not isinstance(match_properties, (tuple, list)):
            match_properties = [match_properties]

        cp.check_list_continuous_properties(match_properties, check_set_value=False)

        # copy and change values to match this stimulus
        match_props = []
        for m in match_properties:
            m = copy(m)
            m.set_value(self._da)
            match_props.append(m)

        da = self._da.copy()
        da.max_array_radius += (extra_space // 2)
        da_sequence = [da]
        error = None

        while (da.prop_numerosity > min_numerosity):
            da = da.number_deviant(change_numerosity=-1)

            if len(match_props)>0:
                da.match(match_props, center_array=True) # TODO center array OK?

            cnt = 0
            while True:
                cnt += 1

                ok, mesg = da.realign()
                if ok:
                    break

                if cnt > 10:
                    error = u"ERROR: realign, " + str(cnt) + ", " + str(da.prop_numerosity)

            da_sequence.append(da)
            if error is not None:
                break

        rtn = DASequence()
        rtn.append_dot_arrays(list(reversed(da_sequence)))
        rtn.method = match_props
        rtn.error = error

        if not inhibit_logging and self.logger is not None:
            self.logger.log(rtn)

        return rtn



class DASequenceGeneratorProcess(Process):

    def __init__(self, max_dot_array, min_numerosity, match_properties, extra_space,
                 n_trials=3, logger=None):

        """
        property: da_sequence, after processes finished
        Event(): data_available
        """

        super(DASequenceGeneratorProcess, self).__init__()


        self.data_available = Event()
        self._data_queue = Queue()
        self._da_sequence = None

        if not isinstance(match_properties, (tuple, list)):
            match_properties = [match_properties]
        cp.check_list_continuous_properties(match_properties)
        self.match_properties = match_properties
        self.max_dot_array = max_dot_array
        self.min_numerosity = min_numerosity
        self.extra_space = extra_space

        if n_trials<1:
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
        super(DASequenceGeneratorProcess, self).join(timeout)

    def _read_queue(self):
        if self._da_sequence is not None:
            return
        if self.is_alive() or not self._data_queue.empty():
            self._da_sequence = self._data_queue.get()

    def run(self):
        cnt = 0
        da_seq = None
        generator = DASequenceGenerator(max_dot_array=self.max_dot_array,
                                        logger=self.logger)

        while cnt<self._n_tryouts:
            cnt += 1
            da_seq = generator.make(match_properties=self.match_properties,
                                    extra_space=self.extra_space,
                                    min_numerosity=self.min_numerosity)
            if da_seq.error is None:
                break
            # print("remix")

        self.data_available.set()
        self._data_queue.put(da_seq)