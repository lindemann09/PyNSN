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


    maxnumber = 200

    logger = GeneratorLogger(log_filename="log/test", override_log_files=True,
                             log_colours=True)
    generator= DotArrayGenerator(
                           max_array_radius= 300,
                           dot_diameter_mean=7,
                           dot_diameter_range=(5, 15),
                           dot_diameter_std=2,
                           dot_colour=(230, 230, 230),
                           minimum_gap=1,
                           logger= logger)

    max_da = generator.make(n_dots=maxnumber, inhibit_logging=True)

    if True:

        x = generator.make(n_dots=100, inhibit_logging=True)
        x.change_colours_random_dots(colours=["lightgreen", "red"],
                                                random_select_ratios=[0.7, 0.3])
        logger.log(x)


    else:
        p = DASequenceGeneratorProcess(max_dot_array=max_da,
                                       match_method=[DASequenceGenerator.CONVEX_HULL,
                                                     DASequenceGenerator.DENSITY_ONLY_AREA],
                                       min_numerosity=10, logger=logger,
                                       extra_space=100)
        p.start()
        p1 = DASequenceGeneratorProcess(max_dot_array=max_da, match_method=DASequenceGenerator.CONVEX_HULL,
                                        min_numerosity=10, logger=logger,
                                       extra_space=100)
        p1.start()
        # x = p.da_sequence


    control.set_develop_mode(True)
    exp = control.initialize()
    control.start(exp, skip_ready_screen=True)


    stim = expyriment_stimulus.ExprimentDotArray(x,
                                          colour_area=(100, 100, 100),
                                           #colour_convex_hull_dots=(255, 0, 0),
                                          #colour_convex_hull_positions=(255, 200, 0),
                                          #colour_center_of_mass=(255, 0, 0),
                                          #  colour_center_of_outer_positions=(0,0,200),
                                                     antialiasing=True)


    stim.present()
    exp.keyboard.wait()



    control.end()