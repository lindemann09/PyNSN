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

    logger = GeneratorLogger(log_filename="log/test", override_log_files=True)
    generator= DotArrayGenerator(
                           field_radius= 300,
                           dot_diameter_mean=7,
                           dot_diameter_range=(5, 15),
                           dot_diameter_std=2,
                           dot_colour=(230, 230, 230),
                           minimum_gap=1,
                           logger= logger)

    max_da = generator.make(n_dots=maxnumber, inhibit_logging=False)
    max_da = generator.make(n_dots=maxnumber, inhibit_logging=False)
    max_da = generator.make(n_dots=maxnumber, inhibit_logging=False)

    if False:

        gen = DASequenceGenerator(max_da, logger=logger)
        x = gen.make(match_methods=DASequenceGenerator.TOTAL_CIRCUMFERENCE, min_numerosity=10)

        #gen = DASequenceGenerator(max_da, logger=logger)
        #y = gen.make(method=DASequenceGenerator.CONVEX_HULL, min_numerosity=10)

    else:
        p = DASequenceGeneratorProcess(max_dot_array=max_da,
                                       match_method=[DASequenceGenerator.CONVEX_HULL,
                                                     DASequenceGenerator.DENSITY_ONLY_AREA],
                                       min_numerosity=10, logger=logger,
                                       sqeeze_factor=0.7)
        p.start()
        p1 = DASequenceGeneratorProcess(max_dot_array=max_da, match_method=DASequenceGenerator.CONVEX_HULL,
                                        min_numerosity=10, logger=logger)
        p1.start()
        p2 = DASequenceGeneratorProcess(max_dot_array=max_da, match_method=DASequenceGenerator.CONVEX_HULL,
                                        min_numerosity=10, logger=logger)
        p2.start()

        p1.join()
        p2.join()
        p.join()
        exit()

    control.set_develop_mode(True)
    exp = control.initialize()
    control.start(exp, skip_ready_screen=True)

    print(x.get_property_string(variable_names=True))

    print(x.dot_arrays[-1].get_property_string(variable_names=True))
    stim = expyriment_stimulus.ExprimentDotArray(x.dot_arrays[-1],
                                          colour_area=(100, 100, 100),
                                           colour_convex_hull_dots=(255, 0, 0),
                                          #colour_convex_hull_positions=(255, 200, 0),
                                          #colour_center_of_mass=(255, 0, 0),
                                          #  colour_center_of_outer_positions=(0,0,200),
                                                     antialiasing=True)


    stim.present()
    exp.keyboard.wait()

    print(x.dot_arrays[0].get_property_string(variable_names=False))
    stim = expyriment_stimulus.ExprimentDotArray(x.dot_arrays[0],
                                          colour_area=(100, 100, 100),
                                           colour_convex_hull_dots=(255, 0, 0),
                                          #colour_convex_hull_positions=(255, 200, 0),
                                          #colour_center_of_mass=(255, 0, 0),
                                          #  colour_center_of_outer_positions=(0,0,200),
                                                     antialiasing=True)


    stim.present()
    exp.keyboard.wait()



    control.end()