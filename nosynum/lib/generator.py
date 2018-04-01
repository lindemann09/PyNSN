from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import os
import sys
import time
import atexit

from . import random_beta, __version__
from Queue import Empty
from multiprocessing import Process, Event, Lock,Queue
from .dot_array import DotArray, DASequence


class DotArrayGenerator(object):

    def __init__(self,
                 stimulus_area_radius,
                 dot_diameter_mean,
                 dot_diameter_range=None,
                 dot_diameter_std=None,
                 dot_picture = None,
                 dot_colour=None,
                 minimum_gap=1,
                 logger=None):

        """Specification of a Random Dot Array

        Parameters:
        -----------
        stimulus_area_radius : int
            the radius of the stimulus area
        n_dots : int
            number of moving dots

        """

        if dot_diameter_std <= 0:
            dot_diameter_std = None
        if dot_diameter_range is not None and \
                (dot_diameter_mean <= dot_diameter_range[0] or
                 dot_diameter_mean >= dot_diameter_range[1] or
                 dot_diameter_range[0] >= dot_diameter_range[1]):
                raise RuntimeError("dot_diameter_mean has to be inside the defined dot_diameter_range")

        self.minimum_gap = minimum_gap
        self.stimulus_area_radius = stimulus_area_radius # TODO rename to array_radius
        self.dot_diameter_range = dot_diameter_range
        self.dot_diameter_mean = dot_diameter_mean
        self.dot_diameter_std = dot_diameter_std
        self.dot_colour = dot_colour
        self.dot_picture = dot_picture
        self.set_logger(logger)


    def set_logger(self, logger):
        self.logger = logger
        if not isinstance(logger, (type(None), GeneratorLogger)):
            raise RuntimeError("logger has to be None or a generator logger") # TODO: check type name


    def make(self, n_dots, inhibit_logging=False):

        rtn = DotArray(max_array_radius = self.stimulus_area_radius,
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
            rtn.append(xy=xy, diameter=diameter,
                        colour=self.dot_colour,
                        picture=self.dot_picture)

        if not inhibit_logging and self.logger is not None:
            self.logger.log(rtn)

        return rtn




class DASequenceGenerator(object):

    MEAN_DIAMETER = "IS"
    CONVEX_HULL = "CH"
    TOTAL_AREA = "TA"
    DENSITY = "DE"
    NO_FITTING = "NF"
    TOTAL_CIRCUMFERENCE = "TC"
    ALL_METHODS = [MEAN_DIAMETER, CONVEX_HULL, TOTAL_AREA, DENSITY, NO_FITTING]

    def __init__(self, max_dot_array, sqeeze_factor=None,  # fixme: squeeze factor needed ?
                 use_convex_hull_positions=False, logger=None):
        """  makes sequence of deviants by subtracting dots

            sqeeze factor: when adapting for convex hull, few point shift excentrically, it is
                  therefore usefull to sqeeze the stimulus before. We do that therefore all stimuli

        Stimulus will be center before making variants

        """
        self._da = []
        self.use_convex_hull_positions = use_convex_hull_positions
        self.set_max_dot_array(max_dot_array=max_dot_array, sqeeze_factor=sqeeze_factor)
        self.set_logger(logger)

    def set_logger(self, logger):
        self.logger = logger
        if not isinstance(logger, (type(None), GeneratorLogger)):
            raise RuntimeError("logger has to be None or a generator logger")  # TODO: check type name

    def set_max_dot_array(self, max_dot_array, sqeeze_factor=None):
        if sqeeze_factor is None:
            sqeeze_factor = 1

        self._da = max_dot_array.copy()
        self._da.fit_convex_hull_area(convex_hull_area=self._da.prop_area_convex_hull_dots * sqeeze_factor)
        self._da.xy -= self._da.center_of_outer_positions # centering

        self._dia = self._da.prop_mean_dot_diameter
        self._cha = self._da.prop_area_convex_hull_dots
        self._dens = self._da.prop_density
        self._total_area = self._da.prop_total_surface_area
        self._circumference = self._da.prop_total_circumference


    def make(self, method, min_numerosity, inhibit_logging=False):
        """Methods takes take , you might use make Process

         returns False is error occured (see self.error)
        """

        da = self._da
        da_sequence = [da]
        error = None

        while (da.prop_numerosity > min_numerosity):
            da = da.number_deviant(change_numerosity=-1)

            if method == DASequenceGenerator.DENSITY:
                da.fit_density(density=self._dens, ratio_area_convex_hull_adaptation=0.5,
                                        use_convex_hull_positions=self.use_convex_hull_positions)
            elif method == DASequenceGenerator.CONVEX_HULL:
                da.fit_convex_hull_area(convex_hull_area=self._cha,
                                        use_convex_hull_positions=self.use_convex_hull_positions)
            elif method == DASequenceGenerator.MEAN_DIAMETER:
                da.fit_mean_dot_diameter(mean_dot_diameter=self._dia)
            elif method == DASequenceGenerator.TOTAL_AREA:
                da.fit_total_area(total_area=self._total_area)
            elif method == DASequenceGenerator.TOTAL_CIRCUMFERENCE:
                da.fit_total_circumference(total_circumference=self._circumference)
            elif method == self.NO_FITTING:
                pass
            else:
                raise Warning("Unknown method {}. Using NO_FITTING.".format(method))

            cnt = 0
            while True:
                cnt += 1
                try:
                    if da.realign():  # Ok if realign not anymore required
                        break
                except:
                    # error = "WARNING: outlier removal, " + str(cnt) + ", " + str(len(da.dots))
                    # print(error)
                    break

                if cnt > 10:
                    error = "ERROR: realign, " + str(cnt) + ", " + str(len(da.dots))

                if error is not None:
                    error = (cnt, error)
                    break

            da_sequence.append(da)
            if error is not None:
                break

        rtn = DASequence()
        rtn.append_dot_arrays(list(reversed(da_sequence)))
        rtn.method = method
        rtn.error = error

        if not inhibit_logging and self.logger is not None:
            self.logger.log(rtn)

        return rtn



class DASequenceGeneratorProcess(Process):

    def __init__(self, max_dot_array, min_numerosity, method,
                 n_trials=3, sqeeze_factor=None, use_convex_hull_positions=False,
                 logger=None):

        """
        property: da_sequence, after processes finished
        Event(): data_available
        """

        super(DASequenceGeneratorProcess, self).__init__()
        self.data_available = Event()
        self._data_queue = Queue()
        self._da_sequence = None
        self.daemon = True

        self.max_dot_array = max_dot_array
        self.method = method
        self.min_numerosity = min_numerosity
        self.sqeeze_factor = sqeeze_factor
        self.use_convex_hull_positions = use_convex_hull_positions

        if n_trials<1:
            self._n_trails = 1
        else:
            self._n_trails = n_trials

        self.logger = logger
        if not isinstance(logger, (type(None), GeneratorLogger)):
            raise RuntimeError("logger has to be None or a generator logger") # TODO: check type name

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
                                        sqeeze_factor=self.sqeeze_factor,
                                        use_convex_hull_positions=self.use_convex_hull_positions,
                                        logger=self.logger)

        while cnt<self._n_trails:
            cnt += 1
            da_seq = generator.make(method=self.method, min_numerosity=self.min_numerosity)
            if da_seq.error is None:
                break
            # print("remix")

        self.data_available.set()
        self._data_queue.put(da_seq)



class GeneratorLogger(Process):

    def __init__(self, log_filename, log_colours=False, num_format="%6.0f",
                 override_log_files=False):
        super(GeneratorLogger, self).__init__()

        self.log_filename_arrays = log_filename + ".array.csv"
        self.log_filename_properties = log_filename + ".prop.csv"
        self.num_format = num_format
        self.log_colours = log_colours
        self.daemon = True

        if override_log_files:
            self.write_mode = "w+"
        else:
            self.write_mode = "a+"

        try:
            os.makedirs(os.path.split(log_filename)[0])
        except: pass

        self.lock = Lock()
        self._quit_event = Event()
        self._varname_written = Event()
        self._new_data_avaiable = Event()
        self._log_queue_array = Queue()
        self._log_queue_prop = Queue()

        self.start()

    def log(self, dot_array_object):

        if isinstance(dot_array_object, (DASequence, DotArray)) and not self._quit_event.is_set():
            txt = dot_array_object.get_property_string(variable_names=not self._varname_written.is_set())
            self._log_queue_prop.put(txt)
            self._new_data_avaiable.set()

            if isinstance(dot_array_object, DotArray):
                txt = dot_array_object.get_csv(hash_column=True, n_dots_column=True,
                                               colour_column=self.log_colours,
                                               num_format=self.num_format,
                                               variable_names=not self._varname_written.is_set())
            else:  # DASequence
                txt = dot_array_object.get_csv(hash_column=True, colour_column=self.log_colours,
                                               num_format=self.num_format,
                                               variable_names=not self._varname_written.is_set())
            self._log_queue_array.put(txt)
            self._new_data_avaiable.set()
            self._varname_written.set()

    def join(self, timeout=None):
        self._quit_event.set()
        self._new_data_avaiable.set()
        super(GeneratorLogger, self).join(timeout)

    def run(self):

        logfile_arrays = open(self.log_filename_arrays, self.write_mode)
        logfile_prop = open(self.log_filename_properties, self.write_mode)

        comment = "# NoSyNum {}, {}, main: {}\n".format(__version__, time.asctime(),
                            os.path.split(sys.argv[0])[1])

        logfile_prop.write(comment)
        logfile_arrays.write(comment)
        while not self._quit_event.is_set():
            self._new_data_avaiable.wait(timeout=1)

            if self._new_data_avaiable.is_set():
                prop_txt = []
                array_txt = []
                # read all queues
                while True:
                    new_data = False
                    try:
                        txt = self._log_queue_prop.get_nowait()
                        prop_txt.append(txt)
                        new_data = True
                    except Empty:
                        pass

                    try:
                        txt = self._log_queue_array.get_nowait()
                        array_txt.append(txt)
                        new_data = True
                    except Empty:
                        pass

                    if not new_data:
                        break

                self._new_data_avaiable.clear()
                # write files
                for txt in prop_txt:
                    logfile_prop.write(txt)
                for txt in array_txt:
                    logfile_arrays.write(txt)

        logfile_arrays.close()
        logfile_prop.close()