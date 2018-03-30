#!/usr/bin/env python
from __future__ import absolute_import, print_function, division

import math
import random
import numpy as np
from time import sleep
from copy import deepcopy
import nosynum
from nosynum import expyriment_stimulus, pil_image
from nosynum.numpy_dot_array import NumpyRandomDotArray
import expyriment
from expyriment import misc, control

cl = misc.Clock()



maxnumber = 200

dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=10,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=2,
                       dot_colour=(230, 230, 230),
                        minium_gap=1)
cl.reset_stopwatch()
max_da = NumpyRandomDotArray(n_dots=maxnumber, dot_array_definition=dot_array_def)
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

control.set_develop_mode(True)
exp = control.initialize()
control.start(exp, skip_ready_screen=True)
pil = pil_image.create(max_da)
expyriment_stimulus.ExprimentPILImage(pil).present()
exp.keyboard.wait()

max_da.realign(minimum_gap=10)
pil = pil_image.create(max_da)
expyriment_stimulus.ExprimentPILImage(pil).present()
exp.keyboard.wait()


control.end()