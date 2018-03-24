#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from builtins import *

import math
import random
import numpy as np
from time import sleep
from copy import deepcopy
import nosynum
from nosynum import expyriment_stimulus, pil_image
from expyriment import misc

cl = misc.Clock()

dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 400,
                       dot_diameter_mean=10,
                       dot_diameter_range=(2, 30),
                       dot_diameter_std=2,
                       dot_colour=(255, 57, 57))

max_da = nosynum.DotArray(n_dots=20, dot_array_definition=dot_array_def)

#print(da.get_array_csv_text(variable_names=True, colour_column=False))

mp = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                        method=nosynum.M_NO_FITTING,
                                        save_images=True)
mp.join()

seq = mp.da_sequence
#print(seq.images)
print(seq.get_csv(variable_names=True, hash_column=True))
print(seq.md5hash)
print(seq.numerosity_idx)
#exp_seq = expyriment_stimulus.ExpyrimentDASequence(seq)
#print("converting")
#cl.reset_stopwatch()
#pil_image.dasequence2images(dot_array_sequence=seq, image_filename="test")
#print(cl.stopwatch_time)

#print(seq.dot_arrays)
#print(seq.get_csv())
#print(seq.numerosity_correlations)




#dal2 = DotArrayList(number_list, subfolder="picts")
#dal2.create_method(method=2, dot_diameter = 10)
#dal.extend(dal2)
#dal2.create_method(method=3, dot_diameter = 10)
#dal.extend(dal2)

#dal.analysis()
#dal.save_arrays_csv_text(filename="dot_arrays.csv")
#dal.save_properties("array_info.csv")

#print "saving"
#dal.save_images(name="test", file_type="JPEG", antialiasing=False)

#import matplotlib.pyplot as plt
#plt.plot(data[:,0], data[:, 4])
#plt.show()
