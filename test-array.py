#!/usr/bin/env python
from __future__ import absolute_import, print_function, division

import math
import random
import numpy as np
from time import sleep
from copy import deepcopy
import nosynum
from nosynum import expyriment_stimulus, DotArrayGenerator, \
    DASequenceGeneratorProcess, DASequenceGenerator, GeneratorLogger
import expyriment
from expyriment import misc, control
import copy

if __name__ == "__main__":
    cl = misc.Clock()

    maxnumber = 20

    logger = GeneratorLogger(log_filename="log/test", override_log_files=True)
    generator= DotArrayGenerator(
                           stimulus_area_radius= 300,
                           dot_diameter_mean=7,
                           dot_diameter_range=(5, 15),
                           dot_diameter_std=2,
                           dot_colour=(230, 230, 230),
                           minimum_gap=1,
                           logger= logger)

    max_da = generator.make(n_dots=maxnumber, inhibit_logging=False)
    max_da = generator.make(n_dots=maxnumber, inhibit_logging=True)

    if False:

        gen = DASequenceGenerator(max_da, logger=logger)
        x = gen.make(method=DASequenceGenerator.CONVEX_HULL, min_numerosity=10)

        gen = DASequenceGenerator(max_da, logger=logger)
        y = gen.make(method=DASequenceGenerator.CONVEX_HULL, min_numerosity=10)

    else:
        p = DASequenceGeneratorProcess(max_dot_array=max_da , method=DASequenceGenerator.CONVEX_HULL,
                                       min_numerosity=10, logger=logger)
        p.start()
        p1 = DASequenceGeneratorProcess(max_dot_array=max_da, method=DASequenceGenerator.CONVEX_HULL,
                                        min_numerosity=10, logger=logger)
        p1.start()
        #print(p.da_sequence.dot_arrays[3].xy)
        #print(p1.da_sequence.dot_arrays[13].xy)
        sleep(0.5)
        x = p.da_sequence
        y = p1.da_sequence

    print(x.md5hash)
    print(y.md5hash)

    exit()
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