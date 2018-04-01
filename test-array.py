#!/usr/bin/env python
from __future__ import absolute_import, print_function, division

import math
import random
import numpy as np
from time import sleep
from copy import deepcopy
import nosynum
from nosynum import expyriment_stimulus, DotArrayGenerator
import expyriment
from expyriment import misc, control

cl = misc.Clock()

maxnumber = 200

generator= DotArrayGenerator(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=7,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=2,
                       dot_colour=(230, 230, 230),
                       minimum_gap=1,
                        log_file=nosynum.GeneratorLogFile(log_filename="log/hhassg.as")
                        )

control.set_develop_mode(True)
exp = control.initialize()
control.start(exp, skip_ready_screen=True)


cl.reset_stopwatch()
max_da = generator.make(n_dots=maxnumber)
stim = expyriment_stimulus.ExprimentDotArray(max_da,
                                      colour_area=(100, 100, 100),
                                       colour_convex_hull_dots=(255, 0, 0),
                                      colour_convex_hull_positions=(255, 200, 0),
                                      colour_center_of_mass=(255, 0, 0),
                                        colour_center_of_outer_positions=(0,0,200),
                                                 antialiasing=True)
stim.preload()
print(cl.stopwatch_time)

stim.present()
exp.keyboard.wait()




control.end()