__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import os
import sys
import time
import gzip
import numpy as np
import atexit

from pynsn import __version__
from pynsn._lib import _misc, _geometry
from pynsn.dot_array_sequence import DASequence
from pynsn._lib._dot_array import DotArray
from pynsn._lib._item_attributes import ItemAttributes

# FIXME REMOVE FILE LOGGING

class LogFile(object):

    def __init__(self, log_filename,
                 log_colours=False,
                 override_log_files=False,
                 zipped=False):

        log_filename_arrays = log_filename + ".array.csv"
        log_filename_properties = log_filename + ".prop.csv"
        self.log_colours = log_colours

        try:
            os.makedirs(os.path.split(log_filename)[0])
        except:
            pass

        if override_log_files:
            write_mode = "wb"
        else:
            write_mode = "ab"

        if zipped:
            self.logfile_arrays = gzip.open(log_filename_arrays + ".gz", write_mode)
            self.logfile_prop = gzip.open(log_filename_properties+ ".gz", write_mode)
        else:
            self.logfile_arrays = open(log_filename_arrays, write_mode)
            self.logfile_prop = open(log_filename_properties, write_mode)

        header = u"# PyNSN {}, {}, main: {}\n".format(__version__, time.asctime(),
                                                      os.path.split(sys.argv[0])[1])
        self.logtext_prop = header
        self.logtext_arrays = header
        self._varname_written = False
        atexit.register(self.close)

    @staticmethod
    def _logging_txt(dot_array_object,
                     variable_names,
                     log_colours):
        """helper function: returns log for dot array and log for properties"""

        if isinstance(dot_array_object, (DASequence, DotArray)):
            da_log = dot_array_object.get_csv(colour_column=log_colours,
                                              variable_names=variable_names)
            is_sequence = isinstance(dot_array_object, DASequence)
            if not is_sequence:
                feat = dot_array_object.features.get_features_dict()
            else:  # DASequence
                feat = dot_array_object.get_features_dict()
            feat_log = _misc.dict_to_csv(feat, variable_names=variable_names,
                                          dict_of_lists=is_sequence)

        else:
            da_log = ""
            feat_log = ""

        return da_log, feat_log

    def log(self, dot_array_object):

        da_log, feat_log = LogFile._logging_txt(dot_array_object=dot_array_object,
                                                variable_names=not self._varname_written,
                                                log_colours=self.log_colours)
        self._varname_written = True
        self.logtext_prop += feat_log
        self.logtext_arrays += da_log

    def __del__(self):
        self.close()

    def close(self):
        if not self.logfile_arrays.closed:
            self.save()
            self.logfile_arrays.close()
            self.logfile_prop.close()

    def save(self):
        self.logfile_arrays.write(self.logtext_arrays.encode("utf-8"))
        self.logfile_prop.write(self.logtext_prop.encode("utf-8"))
        self.logtext_prop = ""
        self.logtext_arrays = ""


class LogFileReader(object): # FIXME really required. Better using files for
    # saving xy

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
        self.hashes = []
        self.num_ids = []
        self.unique_hashes = []

    def load(self):

        with open(self.filename, "r") as fl:
            self.unique_hashes = []
            xy = []
            diameters = []
            colours = []  # TODO: features
            pictures = []
            hashes = []
            num_ids = []
            first_line = True
            for l in fl:
                if l.strip().startswith(self.comment):
                    continue
                if first_line and self.first_line_variable_names:
                    first_line = False
                    continue

                arr = l.split(",")

                hashes.append(arr[0].strip())
                if hashes[-1] not in self.unique_hashes:
                    self.unique_hashes.append(hashes[-1])

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
        self.hashes = np.array(hashes)
        self.num_ids = np.array(num_ids, dtype=np.int)

    def _check_loaded(self):

        if len(self.unique_hashes) == 0:
            raise RuntimeError("Operation not possible. Please first load log file (LogFileReader.load()).")

    def get_object_type(self, hash):

        l = len(self.get_unique_num_ids(hash))
        if l > 1:  # several num_ids or just one
            return DASequence
        elif l == 1:
            return DotArray
        else:
            return None

    def get_unique_num_ids(self, hash):

        self._check_loaded()
        if hash in self.unique_hashes:
            ids = self.hashes == hash
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

        rtn = DotArray(target_array_radius=target_array_radius,
                       minimum_gap=2)

        for x in zip(xy, dia, col, pict):
            rtn.append(xy=x[0], item_diameters=x[1],
                       attributes=ItemAttributes(colour=x[2], picture=x[3]))
        if target_array_radius is None:
            # adjust max_radius if not defined
            radii = _geometry.cartesian2polar(rtn._xy, radii_only=True) + rtn._diameters / 2
            rtn.target_array_radius = int(np.ceil(np.max(radii)))
        return rtn

    def get_object(self, hash, target_array_radius=None):  # todo: save max radius?
        """if target_array_radius=None find minimal required radius

        :returns dot array or dot array sequence"""

        self._check_loaded()
        if hash in self.unique_hashes:

            o_idx = self.hashes == hash
            ot = self.get_object_type(hash)
            if ot == DotArray:
                return self._get_dot_array(o_idx, target_array_radius=target_array_radius)

            elif ot == DASequence:
                rtn = DASequence()

                array_radii = []
                for num_id in self.get_unique_num_ids(hash):
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
