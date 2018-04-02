from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from . import random_beta
from multiprocessing import Process, Event, Queue
from .dot_array import DotArray
from .dot_array_sequence import DASequence
from .files import GeneratorLogger


class DotArrayGenerator(object):

    def __init__(self,
                 max_array_radius,
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
        self.max_array_radius = max_array_radius
        self.dot_diameter_range = dot_diameter_range
        self.dot_diameter_mean = dot_diameter_mean
        self.dot_diameter_std = dot_diameter_std
        self.dot_colour = dot_colour
        self.dot_picture = dot_picture
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
            rtn.append_numpy(xy=xy, diameter=diameter,
                             colour=self.dot_colour,
                             picture=self.dot_picture)

        if not inhibit_logging and self.logger is not None:
            self.logger.log(rtn)

        return rtn




class DASequenceGenerator(object):

    NO_FITTING = "NF"
    CONVEX_HULL = "CH"
    DENSITY = "DE"
    DENSITY_ONLY_CONVEX_HULL ="DE_C"
    DENSITY_ONLY_AREA ="DE_A"
    MEAN_DIAMETER = "MD"
    TOTAL_AREA = "TA"
    TOTAL_CIRCUMFERENCE = "TC"

    ALL_METHODS = [MEAN_DIAMETER, CONVEX_HULL, TOTAL_AREA, DENSITY, NO_FITTING,
                   DENSITY_ONLY_CONVEX_HULL, DENSITY_ONLY_AREA]

    _DEPENDENCIES =[ [CONVEX_HULL, DENSITY, DENSITY_ONLY_CONVEX_HULL],
                     [DENSITY, DENSITY_ONLY_AREA, MEAN_DIAMETER, TOTAL_CIRCUMFERENCE, TOTAL_AREA ],
                     [DENSITY, DENSITY_ONLY_AREA, DENSITY_ONLY_CONVEX_HULL]]

    def __init__(self, max_dot_array, logger=None):
        """  makes sequence of deviants by subtracting dots

            sqeeze factor: when adapting for convex hull, few point shift excentrically, it is
                  therefore usefull to sqeeze the stimulus before. We do that therefore all stimuli

        Stimulus will be center before making variants

        """
        self._da = []
        self.set_max_dot_array(max_dot_array=max_dot_array)
        self.set_logger(logger)

    def set_logger(self, logger):
        self.logger = logger
        if not isinstance(logger, (type(None), GeneratorLogger)):
            raise RuntimeError("logger has to be None or a GeneratorLogger")

    def set_max_dot_array(self, max_dot_array):

        self._da = max_dot_array.copy()
        # auxiliary variables
        self._dia = self._da.prop_mean_dot_diameter
        self._cha = self._da.prop_convex_hull_area
        self._dens = self._da.prop_density
        self._total_area = self._da.prop_total_surface_area
        self._circumference = self._da.prop_total_circumference

    @staticmethod
    def check_match_method_compatibility(match_properties):
        # check compatible combinations
        for dep in DASequenceGenerator._DEPENDENCIES:
            if sum(list(map(lambda x: x in dep, match_properties))) > 1:
                raise RuntimeError("Incompatible properties to match: {}".format(match_properties))
                return False
        return True




    def make(self, match_methods, min_numerosity,
             extra_space, # fitting convex hull and density might result in enlarged arrays
             inhibit_logging=False):
        """Methods takes take , you might use make Process

         returns False is error occured (see self.error)
        """

        if isinstance(match_methods, (tuple, list)):
            DASequenceGenerator.check_match_method_compatibility(match_methods)
        else:
            match_methods = [match_methods]

        da = self._da.copy()
        da.max_array_radius += (extra_space // 2)
        da_sequence = [da]
        error = None

        while (da.prop_numerosity > min_numerosity):
            da = da.number_deviant(change_numerosity=-1)

            for mp in match_methods:
                if mp == DASequenceGenerator.DENSITY:
                    da.match_density(density=self._dens, ratio_convex_hull2area_adaptation=0.5)

                elif mp == DASequenceGenerator.DENSITY_ONLY_AREA:
                    da.match_density(density=self._dens, ratio_convex_hull2area_adaptation=0)

                elif mp == DASequenceGenerator.DENSITY_ONLY_CONVEX_HULL:
                    da.match_density(density=self._dens, ratio_convex_hull2area_adaptation=1)

                elif mp == DASequenceGenerator.CONVEX_HULL:
                    da.match_convex_hull_area(convex_hull_area=self._cha)

                elif mp == DASequenceGenerator.MEAN_DIAMETER:
                    da.match_mean_dot_diameter(mean_dot_diameter=self._dia)

                elif mp == DASequenceGenerator.TOTAL_AREA:
                    da.match_total_surface_area(surface_area=self._total_area)

                elif mp == DASequenceGenerator.TOTAL_CIRCUMFERENCE:
                    da.match_total_circumference(total_circumference=self._circumference)

                elif mp == self.NO_FITTING:
                    pass
                else:
                    raise Warning("Unknown method {}. Using NO_FITTING.".format(mp))

            cnt = 0
            while True:
                cnt += 1

                ok, mesg = da.realign()
                if ok:
                    break

                if cnt > 10:
                    error = "ERROR: realign, " + str(cnt) + ", " + str(da.prop_numerosity)

            da_sequence.append(da)
            if error is not None:
                break

        rtn = DASequence()
        rtn.append_dot_arrays(list(reversed(da_sequence)))
        rtn.method = match_methods
        rtn.error = error

        if not inhibit_logging and self.logger is not None:
            self.logger.log(rtn)

        return rtn



class DASequenceGeneratorProcess(Process):

    def __init__(self, max_dot_array, min_numerosity, match_method, extra_space,
                 n_trials=3, logger=None):

        """
        property: da_sequence, after processes finished
        Event(): data_available
        """

        super(DASequenceGeneratorProcess, self).__init__()


        self.data_available = Event()
        self._data_queue = Queue()
        self._da_sequence = None

        if isinstance(match_method, (tuple, list)):
            DASequenceGenerator.check_match_method_compatibility(match_method)
        self.match_method = match_method
        self.max_dot_array = max_dot_array
        self.min_numerosity = min_numerosity
        self.extra_space = extra_space

        if n_trials<1:
            self._n_trails = 1
        else:
            self._n_trails = n_trials

        self.logger = logger
        if not isinstance(logger, (type(None), GeneratorLogger)):
            raise RuntimeError("logger has to be None or a GeneratorLogger")

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

        while cnt<self._n_trails:
            cnt += 1
            da_seq = generator.make(match_methods=self.match_method,
                                    extra_space=self.extra_space,
                                    min_numerosity=self.min_numerosity)
            if da_seq.error is None:
                break
            # print("remix")

        self.data_available.set()
        self._data_queue.put(da_seq)