#!/usr/bin/env python
from __future__ import absolute_import, print_function, division

import math
import random
import numpy as np
from time import sleep
from copy import deepcopy
import nosynum
from nosynum import expyriment_stimulus, RandomDotArrayGenerator
import expyriment
from expyriment import misc, control

cl = misc.Clock()



maxnumber = 100

generator= RandomDotArrayGenerator(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=7,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=2,
                       dot_colour=(230, 230, 230),
                        minimum_gap=1)
cl.reset_stopwatch()
max_da = generator.make(n_dots=maxnumber)
print(cl.stopwatch_time)
#print(da.get_array_csv_text(variable_names=True, colour_column=False))

#mp = nosynum.MakeDASequenceProcess(max_dot_array=max_da,
#                                        method=nosynum.M_TOTAL_CIRCUMFERENCE,
#                                        min_numerosity=10)
#mp.start()
#seq = mp.da_sequence
#print(seq.get_property_string())
#print(seq.numerosity_correlations)
#print(seq.md5hash)

#print(max_da.get_csv())print(max_da.prop_area_convex_hull_positions)
print(max_da.prop_area_convex_hull_dots)

#exit()

control.set_develop_mode(True)
exp = control.initialize()
control.start(exp, skip_ready_screen=True)

cl.reset_stopwatch()

stim = expyriment_stimulus.ExprimentDotArray(max_da,
                                      colour_area=(100, 100, 100), colour_convex_hull_dots=(255, 0, 0),
                                      colour_convex_hull_positions=(255, 200, 0),
                                       colour_center_of_mass=(255, 0, 0))
print(cl.stopwatch_time)
cl.reset_stopwatch()
stim.preload()
print(cl.stopwatch_time)

stim.present()
exp.keyboard.wait()




control.end()