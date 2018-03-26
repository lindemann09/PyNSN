#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from multiprocessing import Process, Queue, Event
from . import DotArray, DASequence
from .dot_array_sequences import M_NO_FITTING


class TemplateDASequenceProcess(Process): # abstract class

    def __init__(self):
        super(TemplateDASequenceProcess, self).__init__()

        self.sequence_available = Event()
        self._data_queue = Queue()
        self._da_sequence = None

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

    def __init__(self, max_dot_array, method=M_NO_FITTING, n_trials=3,
                 auto_start_process=True):

        """
        property: da_sequence, after processes finished
        Event(): sequence_available
        """

        super(MakeDASequenceProcess, self).__init__()
        self.max_dot_array = max_dot_array
        self.method = method
        self.daemon = True

        if n_trials<1:
            self._n_trails = 1
        else:
            self._n_trails = n_trials

        if isinstance(max_dot_array, DotArray) and auto_start_process:
            self.start()

    def run(self):
        cnt = 0
        da_seq = DASequence()

        while cnt<self._n_trails:
            cnt += 1
            if da_seq.make_by_incrementing(max_dot_array=self.max_dot_array, method=self.method):
                break
            else:
                print("remix")

        if cnt>=self._n_trails:
            raise Warning("Could not fine a solution!")

        self.sequence_available.set()
        self._data_queue.put(da_seq)

