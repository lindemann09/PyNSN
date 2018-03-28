#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from time import sleep
from multiprocessing import Process, Queue, Event
from . import DotArray, DASequence
from .dot_array_sequences import M_NO_FITTING


class TemplateDASequenceProcess(Process): # abstract class

    def __init__(self):
        super(TemplateDASequenceProcess, self).__init__()

        self.data_available = Event()
        self._data_queue = Queue()
        self._da_sequence = None
        self.daemon = True


    @property
    def da_sequence(self):
        self._read_queue()
        return self._da_sequence

    def join(self, timeout=None):
        self._read_queue()
        super(TemplateDASequenceProcess, self).join(timeout)

    def _read_queue(self):

        while self.is_alive():
            try:
                x = self._data_queue.get(timeout=0.1)
            except:
                continue

            self._da_sequence = x
            break


class MakeDASequenceProcess(TemplateDASequenceProcess):

    def __init__(self, max_dot_array, min_numerosity, method=M_NO_FITTING, n_trials=3,
                            sqeeze_factor=None):

        """
        property: da_sequence, after processes finished
        Event(): data_available
        """

        super(MakeDASequenceProcess, self).__init__()
        self.max_dot_array = max_dot_array
        self.method = method
        self.min_numerosity = min_numerosity
        self.sqeeze_factor = sqeeze_factor

        if n_trials<1:
            self._n_trails = 1
        else:
            self._n_trails = n_trials

    def run(self):
        cnt = 0
        da_seq = DASequence()

        while cnt<self._n_trails:
            cnt += 1
            if da_seq.make_by_incrementing(max_dot_array=self.max_dot_array,
                                           method=self.method,
                                           sqeeze_factor=self.sqeeze_factor,
                                           min_numerosity=self.min_numerosity):
                break
            print("remix")

        self.data_available.set()
        self._data_queue.put(da_seq)


class ProcessContainer(object):

    def __init__(self, forerun):
        """forerun: number of processes started in advance,
        if one process is retrieved next process will be started"""
        self.container = []
        self.forerun = forerun

    def add_process(self, process):
        self.container.append(process)
        self._start_processes()

    def _start_processes(self):
        """start process if required"""
        for p in self.container[:self.forerun]:
            if not p.is_alive():
                p.start()

    def pop_processes(self):
        """returns next finished sequence process"""

        pro = self.container.pop(0)
        self._start_processes()
        return pro

    @property
    def running_processes(self):
        return len(list(filter(lambda x:x.is_alive(), self.container)))

    @property
    def data_avaiable(self):
        return len(list(filter(lambda x:x.data_available.is_set(), self.container)))
