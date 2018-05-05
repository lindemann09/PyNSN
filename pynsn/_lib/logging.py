from __future__ import print_function, division, unicode_literals

from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import os
import sys
import time
import numpy as np
from .. import __version__
from . import misc
from .item_attributes import ItemAttributeList
from .dot_array import DotArray
from .dot_array_sequence import DASequence
import atexit


class LogFile(object):

    def __init__(self, log_filename,
                 log_colours=False,
                 num_format="%6.0f",
                 override_log_files=False,
                 properties_different_colour=False):

        self.log_filename_arrays = log_filename + ".array.csv"
        self.log_filename_properties = log_filename + ".prop.csv"
        self.num_format = num_format
        self.log_colours = log_colours
        self.properties_different_colour = properties_different_colour

        if override_log_files:
            write_mode = "w+"
        else:
            write_mode = "a+"

        try:
            os.makedirs(os.path.split(log_filename)[0])
        except:
            pass

        header = u"# PyNSN {}, {}, main: {}\n".format(__version__, time.asctime(),
                                                      os.path.split(sys.argv[0])[1])

        with open(self.log_filename_arrays, write_mode) as logfile_arrays:
            logfile_arrays.write(header)
        with open(self.log_filename_properties, write_mode) as logfile_prop:
            logfile_prop.write(header)
        self.logtext_prop = ""
        self.logtext_arrays = ""

        self._varname_written = False
        atexit.register(self.save)

    @staticmethod
    def _logging_txt(dot_array_object,
                     variable_names,
                     num_format,
                     properties_different_colour,
                     log_colours):
        """helper function: returns log for dot arry and log for properties"""

        if isinstance(dot_array_object, (DASequence, DotArray)):
            is_sequence = isinstance(dot_array_object, DASequence)
            feat = dot_array_object.get_features_dict()
            feat_log = misc.dict_to_csv(feat, variable_names=variable_names,
                                        dict_of_lists=is_sequence)
            if not is_sequence:
                if properties_different_colour:
                    feat = dot_array_object.get_features_split_by_colours()  # todo: check logging different colours
                    if feat is not None:
                        feat_log += feat.get_csv(feat, dict_of_lists=True, variable_names=False)

                da_log = dot_array_object.get_csv(colour_column=log_colours,
                                                  num_format=num_format, variable_names=variable_names)
            else:  # DASequence
                da_log = dot_array_object.get_csv(colour_column=log_colours,
                                                  num_format=num_format, variable_names=variable_names)
        else:
            da_log = ""
            feat_log = ""

        return da_log, feat_log

    def log(self, dot_array_object):

        da_log, feat_log = LogFile._logging_txt(dot_array_object=dot_array_object,
                                                variable_names=not self._varname_written,
                                                num_format=self.num_format,
                                                properties_different_colour=self.properties_different_colour,
                                                log_colours=self.log_colours)
        self._varname_written = True
        self.logtext_prop += feat_log
        self.logtext_arrays += da_log

    def save(self):

        with open(self.log_filename_arrays, "a+") as logfile_arrays:
            logfile_arrays.write(self.logtext_arrays)
        with open(self.log_filename_properties, "a+") as logfile_prop:
            logfile_prop.write(self.logtext_prop)
        self.logtext_prop = ""
        self.logtext_arrays = ""


class LogFileReader(object):
    def __init__(self, filename, colours=False, pictures=False,
                 comment="#", first_line_variable_names=True, zipped=False):  # todo: zip
        self.filename = filename
        self.has_colours = colours
        self.has_pictures = pictures
        self.zipped = zipped
        self.comment = comment
        self.first_line_variable_names = first_line_variable_names
        self.unload()

    def unload(self):
        self.xy = []
        self.diameters = []
        self.colours = []
        self.pictures = []
        self.object_ids = []
        self.num_ids = []
        self.unique_object_ids = []

    def load(self):

        with open(self.filename, "r") as fl:
            self.unique_object_ids = []
            xy = []
            diameters = []
            colours = []  # TODO: features
            pictures = []
            object_ids = []
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
                    colours.append(arr[i].strip())
                    i += 1

                if self.has_pictures:
                    pictures.append(arr[i].strip())

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

    def _get_dot_array(self, idx, target_array_radius):
        """Please do not use this method and use get_object()

        helper method with plausibility check
        """

        xy = self.xy[idx, :]
        dia = self.diameters[idx]
        if self.has_colours:
            col = self.colours[idx]
        else:
            col = [None] * len(xy)
        if self.pictures:
            pict = self.pictures[idx]
        else:
            pict = [None] * len(xy)

        rtn = DotArray(target_array_radius=target_array_radius)
        for x in zip(xy, dia, col, pict):
            rtn.append(xy=x[0], item_diameters=x[1],
                       attributes=ItemAttributeList(colours=x[2], pictures=x[3]))
        if target_array_radius is None:
            # adjust max_radius if not defined
            radii = rtn._cartesian2polar(rtn._xy, radii_only=True) + rtn._diameters / 2
            rtn.target_array_radius = int(np.ceil(np.max(radii)))
        return rtn

    def get_object(self, object_id, target_array_radius=None):  # todo: save max radius?
        """if target_array_radius=None find minimal required radius

        :returns dot array or dot array sequence"""

        self._check_loaded()
        if object_id in self.unique_object_ids:

            o_idx = self.object_ids == object_id
            ot = self.get_object_type(object_id)
            if ot == DotArray:
                return self._get_dot_array(o_idx, target_array_radius=target_array_radius)

            elif ot == DASequence:
                rtn = DASequence()

                array_radii = []
                for num_id in self.get_unique_num_ids(object_id):
                    tmp = o_idx & (self.num_ids == num_id)
                    da = self._get_dot_array(tmp, target_array_radius=target_array_radius)
                    rtn.append_dot_arrays(da)
                    if target_array_radius is None:
                        array_radii.append(da.target_array_radius)

                if target_array_radius is None:
                    # just all area radii to largest
                    tmp = np.max(array_radii)
                    for i in range(len(rtn.dot_arrays)):
                        rtn.dot_arrays[i].target_array_radius = tmp

                return rtn

        return None
