"""
Dot Array Sequence
"""
from __future__ import print_function, division, unicode_literals
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from hashlib import md5
import numpy as np
from .simple_dot_array import DotArrayProperties
from .dot_array import DotArray

class DASequence(object):

    def __init__(self):
        """ docu the use of numerosity_idx see get_array_numerosity
        dot array, please use append_dot_array and delete_array to modify the list

        """

        self.dot_arrays = []
        self.method = None
        self.error = None
        self.numerosity_idx = {}

    def append_dot_arrays(self, arr):
        if isinstance(arr, DotArray):
            arr = [arr]
        self.dot_arrays.extend(arr)
        self.numerosity_idx = {da.prop_numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def delete_dot_arrays(self, array_id):
        self.dot_arrays.pop(array_id)
        self.numerosity_idx = {da.prop_numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def get_array_numerosity(self, number_of_dots):
        """returns array with a particular numerosity"""

        try:
            return self.dot_arrays[self.numerosity_idx[number_of_dots]]
        except:
            return None

    @property
    def min_max_numerosity(self):
        return (self.dot_arrays[0].prop_numerosity, self.dot_arrays[-1].prop_numerosity)

    @property
    def object_id(self):
        """meta hash of all of csv (n_dots, counter, position, diameter only)"""

        m = md5()
        for da in self.dot_arrays:
            m.update(da.object_id.encode("UTF-8"))
        return m.hexdigest()[:DotArray.OBJECT_ID_LENGTH]

    def get_properties(self):
        """named tuple with arrays"""
        rtn = DotArrayProperties._make_arrays()

        for da in self.dot_arrays:
            prop = da.get_properties()
            for i in range(len(rtn)):
                rtn[i].append(prop[i])
        return rtn

    def get_numerosity_correlations(self):
        prop = self.get_properties()
        cor = np.corrcoef(np.round(prop.np_array, 2), rowvar=False)
        cor = cor[0, :]
        names = prop.property_names
        rtn = {}
        for x in range(1, len(cor)):
            rtn[names[x]] = cor[x]
        return rtn

    def __str__(self):
        return self.get_csv()

    def get_csv(self, num_format="%7.2f", variable_names=True, colour_column=False,
                picture_column=False, object_id_column=True):

        rtn = ""
        tmp_var_names = variable_names
        for da in self.dot_arrays:
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
