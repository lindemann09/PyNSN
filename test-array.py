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



maxnumber = 200

generator= RandomDotArrayGenerator(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=7,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=2,
                       dot_colour=(230, 230, 230),
                        minimum_gap=1)

control.set_develop_mode(True)
exp = control.initialize()
control.start(exp, skip_ready_screen=True)

max_da = generator.make(n_dots=maxnumber)
print(max_da.prop_area_convex_hull_dots)
stim = expyriment_stimulus.ExprimentDotArray(max_da,
                                      colour_area=(100, 100, 100),
                                       colour_convex_hull_dots=(255, 0, 0),
                                      colour_convex_hull_positions=(255, 200, 0),
                                      colour_center_of_mass=(255, 0, 0))
stim.present()
exp.keyboard.wait()

max_da.minimum_gap = 10
print(max_da.realign())
stim = expyriment_stimulus.ExprimentDotArray(max_da,
                                      colour_area=(100, 100, 100),
                                       colour_convex_hull_dots=(255, 0, 0),
                                      colour_convex_hull_positions=(255, 200, 0),
                                      colour_center_of_mass=(255, 0, 0))
stim.present()
print(max_da.prop_area_convex_hull_dots)

exp.keyboard.wait()


control.end()