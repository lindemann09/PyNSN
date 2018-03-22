#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from builtins import *

import os
import pickle
import gzip
from time import sleep
from multiprocessing import Process, Queue, queues, cpu_count
import numpy as np

from .dot_array import DotArray

M_ITEM_SIZE = 2
M_CONVEX_HULL = 3
M_TOTAL_AREA = 4
M_DENSITY = 5

EXTENSION = ".das"

def make_da_sequence(max_dot_array, method):
    """makes sequence of deviants by subtracting dots

    error queue: multiprocess.Queue for multiprocessing
        if runtime error occured it writes with tuple (dot_array, method)
    """

    da = max_dot_array.copy()
    da_sequence = [da]
    error = False
    cha = da.convex_hull_area
    dens = da.density
    dot_diameter = da.mean_dot_diameter
    total_area = da.total_area
    for x in range(len(da.dots)-10):
        da = da.number_deviant(change_numerosity=-1)
        cnt = 0

        while True:
            cnt += 1

            if method == M_DENSITY:
                da.fit_density(density=dens, ratio_area_convex_hull_adaptation=0.5)
            elif method == M_CONVEX_HULL:
                da.fit_convex_hull_area(convex_hull_area=cha)
            elif method == M_ITEM_SIZE:
                da.fit_mean_item_size(mean_dot_diameter=dot_diameter)
            elif method == M_TOTAL_AREA:
                da.fit_total_area(total_area=total_area)
            else:
                pass

            try:
                if not da.realign():  # Ok if realign not anymore required
                    break
            except:
                error = "ERROR: outlier, " + str(cnt) + ", " + str(len(da.dots))

            if cnt>100:
                error = "ERROR: realign, " + str(cnt) + ", " + str(len(da.dots))

            if error != False:
                print(error)
                error = (cnt, error)
                break

        da_sequence.append(da)
        if error != False:
            break

    return (reversed(da_sequence), error)


class MakeDASequenceProcess(Process):

    def __init__(self, max_dot_array, method):
        super(MakeDASequenceProcess, self).__init__()
        self.max_dot_array = max_dot_array
        self.method = method
        self._error_queue = Queue()
        self._da_sequence_queue = Queue()
        self._error = None
        self._da_sequence = []

    def run(self):
        da_seq, error = make_da_sequence(max_dot_array=self.max_dot_array, method=self.method)
        self._error_queue.put(error)
        for x in da_seq:
            self._da_sequence_queue.put(x)

    def get_results(self):
        self.join()
        if not self._error_queue.empty(): # new data
            self._error = self._error_queue.get_nowait()
            while not self._da_sequence_queue.empty():
                self._da_sequence.append(self._da_sequence_queue.get_nowait())

        return( (self._da_sequence, self._error) )

def make_multiple_dot_array_sequences(n_versions,
                                      methods,
                                      n_dots,
                                      stimulus_area_radius,
                                      dot_diameter_mean,
                                      dot_diameter_range=None,
                                      dot_diameter_std=None,
                                      dot_colour=None,
                                      minium_gap=1,
                                      sqeeze_factor=.9,
                                      subfolder = "stimuli",
                                      extension = EXTENSION,
                                      n_processes = cpu_count(),
                                      override_existing_arrays=False,
                                      save_erroneous_arrays=False):
    """!!!! # FIXME
    returns also feedback

    """

    try:
        os.mkdir(subfolder)
    except:
        pass

    # make processes
    processes = []
    for cnt in range(n_versions):
        for method in methods:
            filename = os.path.join(subfolder, method+"-"+str(cnt)+ extension)
            if override_existing_arrays or not os.path.isfile(filename):
                print("preparing "+filename)
                da = DotArray(n_dots=n_dots, stimulus_area_radius=stimulus_area_radius,
                              dot_diameter_mean=dot_diameter_mean, dot_diameter_range=dot_diameter_range,
                              dot_diameter_std=dot_diameter_std, dot_colour=dot_colour, minium_gap=minium_gap)

                if sqeeze_factor<1 and sqeeze_factor>0:
                    da.fit_convex_hull_area(convex_hull_area=da.convex_hull_area * sqeeze_factor)

                error_queue = Queue()
                sequence_queue = Queue()
                p = Process(target = make_da_sequence, args=(da, method, sequence_queue, error_queue))
                processes.append((p, error_queue, sequence_queue) )

    # run processes
    running=[]
    while len(processes)>0:
        # remove dead processes from running
        for p in running:
            if not p.is_alive():
                running.remove(p)

        # move process to running and start if possible
        if len(running) < n_processes:
            running.append(processes.pop(0))
            running[-1].start()
        else:
            sleep(1)

    # wait running processes are to terminate
    for p in running:
        p.join()

    realign_error_list = []
    while not error_queue.empty():
        realign_error_list.append(error_queue.get_nowait())

    return realign_error_list


class DotArraySequence(object):

    def __init__(self, path):
        with gzip.open(path, 'rb') as fl:
            self.da_sequence = pickle.load(file=fl)

        self._n_dots = list(map(lambda x: len(x.dots), self.da_sequence))

    def get_number(self, number):
        try:
            idx = self._n_dots.index(number)
            return self.da_sequence[idx]
        except:
            return None

    @property
    def property_names(self):
        return self.da_sequence[0].property_names

    @property
    def properties(self):
        rtn = []
        for da in self.da_sequence:
            rtn. append(da.properties)
        return np.array(rtn)

    @property
    def property_string(self):
        rtn = self.property_names + "\n"
        for da in self.da_sequence:
            rtn += str(da.properties).replace("[", "").replace("]", "\n")

        return rtn[:-1]
