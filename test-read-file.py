#!/usr/bin/env python
from __future__ import absolute_import, print_function, division

import math
import random
import numpy as np
from time import sleep
from copy import deepcopy
from scipy.spatial import ConvexHull

import pynsn
from pynsn import expyriment_stimulus, DotArrayGenerator, \
    DASequenceMakeProcess, DASequenceGenerator, GeneratorLogger, DotArray, DASequence, colours
import expyriment
from expyriment import misc, control
import copy


if __name__ == "__main__":
    cl = misc.Clock()

    reader = pynsn.LogFileReader("log/test.array.csv", colours=True)
    reader.load()
    for x in reader.unique_object_ids:
        print((x, reader.get_object_type(x) ))

    x = reader.get_object("64860e3f", max_array_radius=None) #f24e98ab, 73d95768
    reader.unload()
    print(x.get_csv(colour_column=True))

    control.set_develop_mode(True)
    exp = control.initialize()
    control.start(exp, skip_ready_screen=True)

    #print(x.get_property_string(variable_names=True))

    #print(x.dot_arrays[-1].get_property_string(variable_names=True))
    stim = expyriment_stimulus.ExprimentDotArray(x,
                                          colour_area=(100, 100, 100),
                                           colour_convex_hull_dots=(255, 0, 0),
                                          #colour_convex_hull_positions=(255, 200, 0),
                                          #colour_center_of_mass=(255, 0, 0),
                                          #  colour_center_of_outer_positions=(0,0,200),
                                                     antialiasing=True)


    stim.present()
    exp.keyboard.wait()


    control.end()