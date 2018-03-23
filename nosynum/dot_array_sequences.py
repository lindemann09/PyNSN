#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from builtins import *

from multiprocessing import Process, Queue, Event
import numpy as np

from . import sequence_methods

class MakeDASequenceProcess(Process):

    def __init__(self, max_dot_array, method, n_trails=3):
        super(MakeDASequenceProcess, self).__init__()
        self.max_dot_array = max_dot_array
        self.method = method
        self.sequence_available = Event()
        self._data_queue = Queue()
        self._da_sequence = None
        self.daemon = True
        if n_trails<1:
            self._n_trails = 1
        else:
            self._n_trails = n_trails

        self.start()

    @property
    def da_sequence(self):
        self._get_results()
        return self._da_sequence

    def run(self):
        cnt = 0
        da_seq = DASequence()

        while cnt<self._n_trails:
            cnt += 1
            if da_seq.make_by_incrementing(max_dot_array=self.max_dot_array, method=self.method):
                break

        self.sequence_available.set()
        self._data_queue.put(da_seq)

    def join(self, timeout=None):
        self._get_results()
        super(MakeDASequenceProcess, self).join(timeout)

    def _get_results(self):

        while self.is_alive():
            try:
                x = self._data_queue.get(timeout=0.1)
            except:
                continue

            self._da_sequence = x
            break



class DASequence(object):

    def __init__(self):
        self.da_sequence = []
        self.method = None
        self.error = None
        self._n_dots = []

    def get_number(self, number):
        try:
            idx = self._n_dots.index(number)
            return self.da_sequence[idx]
        except:
            return None

    @property
    def max_dot_array(self):
        return self.da_sequence[-1]

    @property
    def property_names(self):
        return self.da_sequence[0].property_names

    @property
    def properties(self):
        rtn = []
        for da in self.da_sequence:
            rtn.append(da.properties)
        return np.array(rtn)

    @property
    def property_string(self):
        rtn = self.property_names + "\n"
        for da in self.da_sequence:
            rtn += str(da.properties).replace("[", "").replace("]", "\n")

        return rtn[:-1]


    def make_by_incrementing(self, max_dot_array, method):
        """makes sequence of deviants by subtracting dots

        returns False is error occured (see self.error)
        """

        da = max_dot_array.copy()

        da_sequence = [da]
        error = None
        cha = da.convex_hull_area
        dens = da.density
        dot_diameter = da.mean_dot_diameter
        total_area = da.total_area
        for x in range(len(da.dots)-10):
            da = da.number_deviant(change_numerosity=-1)
            cnt = 0

            while True:
                cnt += 1

                if method == sequence_methods.DENSITY:
                    da.fit_density(density=dens, ratio_area_convex_hull_adaptation=0.5)
                elif method == sequence_methods.CONVEX_HULL:
                    da.fit_convex_hull_area(convex_hull_area=cha)
                elif method == sequence_methods.ITEM_SIZE:
                    da.fit_mean_item_size(mean_dot_diameter=dot_diameter)
                elif method == sequence_methods.TOTAL_AREA:
                    da.fit_total_area(total_area=total_area)
                elif method == sequence_methods.NO_FITTING:
                    pass
                else:
                    method == sequence_methods.NO_FITTING
                    raise Warning("Unknown method {}. Using NO_FITTING.".format(method))

                try:
                    if not da.realign():  # Ok if realign not anymore required
                        break
                except:
                    error = "ERROR: outlier, " + str(cnt) + ", " + str(len(da.dots))

                if cnt>100:
                    error = "ERROR: realign, " + str(cnt) + ", " + str(len(da.dots))

                if error is not None:
                    error = (cnt, error)
                    break

            da_sequence.append(da)
            if error is not None:
                break

        self.da_sequence = list(reversed(da_sequence))
        self._n_dots = list(map(lambda x: len(x.dots), self.da_sequence))
        self.method = method
        self.error = error

        return error==False