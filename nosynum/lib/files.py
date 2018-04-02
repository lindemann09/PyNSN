from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import os
import sys
import time

from . import __version__
from multiprocessing import Process, Event, Queue
from .dot_array import DotArray
from .dot_array_sequence import DASequence


class GeneratorLogger(Process):

    def __init__(self, log_filename, log_colours=False, num_format="%6.0f",
                 override_log_files=False):
        super(GeneratorLogger, self).__init__()

        self.log_filename_arrays = log_filename + ".array.csv"
        self.log_filename_properties = log_filename + ".prop.csv"
        self.num_format = num_format
        self.log_colours = log_colours

        if override_log_files:
            self.write_mode = "w+"
        else:
            self.write_mode = "a+"

        try:
            os.makedirs(os.path.split(log_filename)[0])
        except: pass

        self._quit_event = Event()
        self._varname_written = Event()
        self._new_data_avaiable = Event()
        self._log_queue_array = Queue()
        self._log_queue_prop = Queue()

        #atexit.register(self.join)
        self.start()

    def log(self, dot_array_object):

        if isinstance(dot_array_object, (DASequence, DotArray)) and not self._quit_event.is_set():
            txt = dot_array_object.get_property_string(variable_names=not self._varname_written.is_set())
            self._log_queue_prop.put(txt)
            self._new_data_avaiable.set()

            if isinstance(dot_array_object, DotArray):
                txt = dot_array_object.get_csv(colour_column=self.log_colours,
                                               num_format=self.num_format,
                                               variable_names=not self._varname_written.is_set())
            else:  # DASequence
                txt = dot_array_object.get_csv(colour_column=self.log_colours,
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
                    except:
                        pass

                    try:
                        txt = self._log_queue_array.get_nowait()
                        array_txt.append(txt)
                        new_data = True
                    except:
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

class LogFileReader(object):

    def __index__(self):
        pass