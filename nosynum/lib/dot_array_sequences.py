#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

import os
import numpy as np

from .dot_array import my_md5_hash

M_ITEM_SIZE = "IS"
M_CONVEX_HULL = "CH"
M_TOTAL_AREA = "TA"
M_DENSITY = "DE"
M_NO_FITTING = "NF"
M_TOTAL_CIRCUMFERENCE = "TC"

ALL_METHODS = [M_ITEM_SIZE, M_CONVEX_HULL, M_TOTAL_AREA, M_DENSITY, M_NO_FITTING]

def is_method(value):
    return value in ALL_METHODS

class DASequence(object):

    def __init__(self):
        self.dot_arrays = []
        self.method = None
        self.error = None
        self.numerosity_idx = {}

    def get_array_numerosity(self, number_of_dots):
        """returns array with a particular numerosity"""

        try:
            return self.dot_arrays[self.numerosity_idx[number_of_dots]]
        except:
            return None

    @property
    def min_max_numerosity(self):
        return (len(self.dot_arrays[0].dots), len(self.dot_arrays[-1].dots))

    @property
    def md5hash(self):
        """md5_hash of csv (n_dots, counter, position, diameter only)"""

        csv = self.get_csv(num_format="%7.2f", hash_column=False,
                           variable_names=False,
                           colour_column=False, picture_column=False)
        return my_md5_hash(csv)

    @property
    def property_names(self):
        return self.dot_arrays[0].property_names

    @property
    def properties(self):
        rtn = []
        for da in self.dot_arrays:
            rtn.append(da.properties)
        return np.array(rtn)

    @property
    def property_correlations(self):
        return np.corrcoef(np.round(self.properties,2), rowvar=False)

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


    def get_property_string(self, varnames=False):
        rtn =""
        if varnames:
            rtn += "hash, " + ", ".join(self.property_names) + "\n"
        hash = self.md5hash
        for da in self.dot_arrays:
            rtn += hash + "," + str(da.properties).replace("[", "").replace("]", "\n")

        return rtn[:-1]

    def get_csv(self, num_format="%7.2f", hash_column=False,
                variable_names=True, colour_column=False, picture_column=False):

        rtn = ""
        tmp_var_names = variable_names
        for da in self.dot_arrays:
            rtn += da.get_csv(n_dots_column=True, hash_column=False,
                              variable_names=tmp_var_names,
                              num_format=num_format, colour_column=colour_column,
                              picture_column=picture_column)
            tmp_var_names = False

        if hash_column:
            hash = self.md5hash
            rtn2 = ""
            tmp_var_names = variable_names
            for l in rtn.split("\n"):
                if tmp_var_names:
                    rtn2 += "hash," + l + "\n"
                    tmp_var_names = False
                elif len(l)>0:
                    rtn2 += "{},{}\n".format(hash, l)
            return rtn2
        else:
            return rtn

    def make_by_incrementing(self, max_dot_array, method, min_numerosity, sqeeze_factor=None):
        """makes sequence of deviants by subtracting dots

        sqeeze factor: when adapting for convex hull, few point shift excentrically, it is
                      therefore need to sqeeze the stimulus before. We do that therefore all stimuli

        returns False is error occured (see self.error)
        """

        if sqeeze_factor is None:
            sqeeze_factor = 1

        da = max_dot_array.copy()
        da.fit_convex_hull_area(convex_hull_area=da.prop_area_convex_hull_positions * sqeeze_factor)
        da_sequence = [da]

        error = None
        cha = da.prop_area_convex_hull_positions
        dens = da.prop_density
        total_area = da.prop_total_area
        circumference = da.prop_total_circumference
        for x in range(len(da.dots)-min_numerosity):
            da = da.number_deviant(change_numerosity=-1)

            if method == M_DENSITY:
                da.fit_density(density=dens, ratio_area_convex_hull_adaptation=0.5)
            elif method == M_CONVEX_HULL:
                da.fit_convex_hull_area(convex_hull_area=cha)
            elif method == M_ITEM_SIZE:
                pass
            elif method == M_TOTAL_AREA:
                da.fit_total_area(total_area=total_area)
            elif method == M_TOTAL_CIRCUMFERENCE:
                da.fit_total_circumference(total_circumference=circumference)
            elif method == M_NO_FITTING:
                pass
            else:
                raise Warning("Unknown method {}. Using NO_FITTING.".format(method))

            cnt = 0
            while True:
                cnt += 1

                try:
                    if not da.realign():  # Ok if realign not anymore required
                        break
                except:
                    #error = "WARNING: outlier removal, " + str(cnt) + ", " + str(len(da.dots))
                    #print(error)
                    break


                if cnt>100:
                    error = "ERROR: realign, " + str(cnt) + ", " + str(len(da.dots))

                if error is not None:
                    error = (cnt, error)
                    break

            da_sequence.append(da)
            if error is not None:
                break

        self.dot_arrays = list(reversed(da_sequence))
        self.numerosity_idx = {len(da.dots):idx for idx, da in enumerate(self.dot_arrays)}
        self.method = method
        self.error = error

        return error is None
