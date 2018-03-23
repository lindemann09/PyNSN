#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from builtins import *

import numpy as np

M_ITEM_SIZE = 2
M_CONVEX_HULL = 3
M_TOTAL_AREA = 4
M_DENSITY = 5
M_NO_FITTING = 6

ALL_METHODS = [M_ITEM_SIZE, M_CONVEX_HULL, M_TOTAL_AREA, M_DENSITY, M_NO_FITTING]

def is_method(value):
    return value in ALL_METHODS

class DASequence(object):

    def __init__(self):
        self.dot_arrays = []
        self.method = None
        self.error = None
        self._n_dots = []

    def get_number(self, number):
        try:
            idx = self._n_dots.index(number)
            return self.dot_arrays[idx]
        except:
            return None

    @property
    def max_dot_array(self):
        return self.dot_arrays[-1]

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
        return np.corrcoef(self.properties, rowvar=False)

    @property
    def numerosity_correlations(self):
        cor = self.property_correlations[0, :]
        rtn = {}
        for x in range(1, len(cor)):
            rtn[self.property_names[x]] = cor[x]
        return rtn


    @property
    def property_string(self):
        rtn = ", ".join(self.property_names) + "\n"
        for da in self.dot_arrays:
            rtn += str(da.properties).replace("[", "").replace("]", "\n")

        return rtn[:-1]

    def get_csv(self, num_format="%7.2f", variable_names=True, colour_column=False, picture_column=False):

        rtn = ""
        for da in self.dot_arrays:
            rtn += da.get_csv(n_dots_column=True, variable_names=variable_names,
                              num_format=num_format, colour_column=colour_column,
                              picture_column=picture_column)
            variable_names = False
        return rtn


    def make_by_incrementing(self, max_dot_array, method):
        """makes sequence of deviants by subtracting dots

        returns False is error occured (see self.error)
        """

        da = max_dot_array.copy()

        da_sequence = [da]
        error = None
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
                elif method == M_NO_FITTING:
                    pass
                else:
                    method == sequence_methods.NO_FITTING
                    raise Warning("Unknown method {}. Using NO_FITTING.".format(method))

                try:
                    if not da.realign():  # Ok if realign not anymore required
                        break
                except:
                    error = "ERROR: outlier, " + str(cnt) + ", " + str(len(da.dots))

                if cnt>100:
                    error = "ERROR: realign, " + str(cnt) + ", " + str(len(da.dots))

                if error is not None:
                    error = (cnt, error)
                    break

            da_sequence.append(da)
            if error is not None:
                break

        self.dot_arrays = list(reversed(da_sequence))
        self._n_dots = list(map(lambda x: len(x.dots), self.dot_arrays))
        self.method = method
        self.error = error

        return error==False