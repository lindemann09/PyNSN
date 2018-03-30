#!/usr/bin/env python
from __future__ import absolute_import, print_function, division

import math
import random
import numpy as np
from time import sleep
from copy import deepcopy
import nosynum
from nosynum import expyriment_pil_image, pil_image, RandomDotArrayCreater
import expyriment
from expyriment import misc, control

cl = misc.Clock()



maxnumber = 100

dot_array_creator= RandomDotArrayCreater(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=7,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=2,
                       dot_colour=(230, 230, 230),
                        minimum_gap=1)
cl.reset_stopwatch()
max_da = dot_array_creator.make(n_dots=maxnumber)
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

#print(max_da.get_csv())
print(max_da.prop_area_convex_hull_positions)
print(max_da.prop_area_convex_hull_dots)

#exit()

control.set_develop_mode(True)
exp = control.initialize()
control.start(exp, skip_ready_screen=True)
pil = pil_image.create(max_da, colour_area=(100, 100, 100), colour_convex_hull_dots=(255, 0, 0),
                       colour_convex_hull_positions=(255, 200, 0),
                       colour_center_of_mass=(255, 0, 0))
expyriment_pil_image.ExprimentPILImage(pil).present()
exp.keyboard.wait()

max_da.realign(minimum_gap=1)
pil = pil_image.create(max_da, colour_area=(100, 100, 100), convex_hull_colour=(255, 0, 0))
expyriment_pil_image.ExprimentPILImage(pil).present()
#exp.keyboard.wait()


control.end()