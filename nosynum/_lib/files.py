from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import os
import sys
import time
import numpy as np

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

    def join(self, timeout=10):
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

    def __init__(self, filename, colours=False, pictures=False,
                 comment = "#", first_line_variable_names=True, zipped=False): #todo: zip
        self.filename = filename
        self.has_colours= colours
        self.has_pictures= pictures
        self.zipped = zipped
        self.comment = comment
        self.first_line_variable_names = first_line_variable_names
        self.unload()

    def unload(self):
        self.xy = []
        self.diameters = []
        self.colours =[]
        self.pictures = []
        self.object_ids = []
        self.num_ids = []
        self.unique_object_ids = []


    def load(self):

        with open(self.filename, "r") as fl:
            self.unique_object_ids = []
            xy = []
            diameters = []
            colours = []
            pictures = []
            object_ids= []
            num_ids = []
            first_line = True
            for l in fl:
                if l.strip().startswith(self.comment):
                    continue
                if first_line and self.first_line_variable_names:
                    first_line = False
                    continue

                arr = l.split(",")

                object_ids.append(arr[0].strip())
                if object_ids[-1] not in self.unique_object_ids:
                    self.unique_object_ids.append(object_ids[-1])

                num_ids.append(int(arr[1]))
                xy.append([float(arr[2]), float(arr[3])])
                diameters.append(float(arr[4]))
                i = 5
                if self.has_colours:
                    c = []
                    for x in range(3):
                        c.append(float(arr[i+x]))
                        i += 1
                    colours.append(c)

                if self.has_pictures:
                    pictures.append(arr[i])

        self.xy = np.array(xy)
        self.diameters = np.array(diameters)
        self.colours = np.array(colours)
        self.pictures = np.array(pictures)
        self.object_ids = np.array(object_ids)
        self.num_ids = np.array(num_ids, dtype=np.int)

    def _check_loaded(self):

        if len(self.unique_object_ids) == 0:
            raise RuntimeError("Operation not possible. Please first load log file (LogFileReader.load()).")


    def get_object_type(self, object_id):

        l = len(self.get_unique_num_ids(object_id))
        if l > 1:  # several num_ids or just one
            return DASequence
        elif l == 1:
            return DotArray
        else:
            return None


    def get_unique_num_ids(self, object_id):

        self._check_loaded()
        if object_id in self.unique_object_ids:
            ids = self.object_ids == object_id
            return np.sort(np.unique(self.num_ids[ids]))
        else:
            return np.array([])


    def _get_dot_array(self, idx, max_array_radius):
        """Please do not use and use get_object()

        helper function with plausibility check
        """

        xy = self.xy[idx, :]
        dia = self.diameters[idx]
        if self.has_colours:
            col = self.colours[idx, :]
        else:
            col = [None] * len(xy)
        if self.pictures:
            pict = self.pictures[idx]
        else:
            pict = [None] * len(xy)

        rtn = DotArray(max_array_radius=max_array_radius)
        for x in zip(xy, dia, col, pict):
            rtn.append_numpy(xy=x[0], diameter=x[1], colour=x[2], picture=x[3])

        return rtn


    def get_object(self, object_id, max_array_radius): #todo: save max radius?

        self._check_loaded()
        if object_id in self.unique_object_ids:

            o_idx = self.object_ids == object_id
            ot = self.get_object_type(object_id)
            if ot == DotArray:
                return self._get_dot_array(o_idx, max_array_radius=max_array_radius)

            elif ot == DASequence:
                rtn = DASequence()

                for num_id in self.get_unique_num_ids(object_id):
                    tmp = o_idx & (self.num_ids == num_id)
                    da = self._get_dot_array(tmp, max_array_radius=max_array_radius)
                    rtn.append_dot_arrays(da)

                return rtn


        return None
