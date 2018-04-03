"""
Dot Array Sequence
"""
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'


import numpy as np
from .dot_array import DotArray, short_md5_hash


class DASequence(object):

    def __init__(self):
        """ docu the use of numerosity_idx see get_array_numerosity"""

        self._dot_arrays = []
        self.method = None
        self.error = None
        self.numerosity_idx = {}

    @property
    def dot_arrays(self):
        """ dot array, please use append_dot_array and delete_array to modify the list"""
        return self._dot_arrays

    def append_dot_arrays(self, arr):
        if isinstance(arr, DotArray):
            arr = [arr]
        self._dot_arrays.extend(arr)
        self.numerosity_idx = {da.prop_numerosity: idx for idx, da in enumerate(self._dot_arrays)}

    def delete_dot_arrays(self, array_id):
        self._dot_arrays.pop(array_id)
        self.numerosity_idx = {da.prop_numerosity: idx for idx, da in enumerate(self._dot_arrays)}

    def get_array_numerosity(self, number_of_dots):
        """returns array with a particular numerosity"""

        try:
            return self._dot_arrays[self.numerosity_idx[number_of_dots]]
        except:
            return None

    @property
    def min_max_numerosity(self):
        return (self._dot_arrays[0].prop_numerosity, self._dot_arrays[-1].prop_numerosity)

    @property
    def object_id(self):
        """md5_hash of csv (n_dots, counter, position, diameter only)"""

        csv = self.get_csv(num_format="%7.2f", object_id_column=False,
                           variable_names=False,
                           colour_column=False, picture_column=False)
        return short_md5_hash(csv, DotArray.OBJECT_ID_LENGTH)

    @property
    def property_names(self):
        return self._dot_arrays[0].property_names

    @property
    def properties(self):
        rtn = []
        for da in self._dot_arrays:
            rtn.append(da.properties)
        return np.array(rtn)

    @property
    def property_correlations(self):
        return np.corrcoef(np.round(self.properties, 2), rowvar=False)

    @property
    def variances(self):
        return np.var(self.properties, axis=0)

    @property
    def numerosity_correlations(self):
        cor = self.property_correlations[0, :]
        rtn = {}
        for x in range(1, len(cor)):
            rtn[self.property_names[x]] = cor[x]
        return rtn

    def get_property_string(self, variable_names=True):
        rtn = ""
        if variable_names:
            rtn += "object_id, " + ", ".join(self.property_names) + "\n"
        obj_id = self.object_id
        for da in self._dot_arrays:
            rtn += obj_id + "," + str(da.properties).replace("[", "").replace("]", "\n")

        return rtn

    def __str__(self):
        return self.get_csv()

    def get_csv(self, num_format="%7.2f", variable_names=True, colour_column=False,
                picture_column=False, object_id_column=True):

        rtn = ""
        tmp_var_names = variable_names
        for da in self._dot_arrays:
            rtn += da.get_csv(num_idx_column=True, object_id_column=False,
                              variable_names=tmp_var_names,
                              num_format=num_format, colour_column=colour_column,
                              picture_column=picture_column)
            tmp_var_names = False

        if object_id_column:
            obj_id = self.object_id
            rtn2 = ""
            tmp_var_names = variable_names
            for l in rtn.split("\n"):
                if tmp_var_names:
                    rtn2 += "object_id," + l + "\n"
                    tmp_var_names = False
                elif len(l) > 0:
                    rtn2 += "{},{}\n".format(obj_id, l)
            return rtn2
        else:
            return rtn
