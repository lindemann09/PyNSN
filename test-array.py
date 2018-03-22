#!/usr/bin/env python

import math
import random
import numpy as np
import nosynum
from expyriment import misc
c = misc.Clock()

dot_diameter_range=(2, 20)
number_list = range(10, 200)

print "creating "
da = nosynum.DotArray(n_dots=100, stimulus_area_radius= 400,
                       dot_diameter_mean=10,
                       dot_diameter_range=(2, 30),
                       dot_diameter_std=2)
print(da.get_array_csv_text())


x = nosynum.make_da_sequence(max_dot_array=da, method=nosynum.M_CONVEX_HULL)
print(x)
#p = nosynum.MakeDASequenceProcess(max_dot_array=da, method=nosynum.M_CONVEX_HULL)
p.start()
x, error = p.get_results()
print(error)
print(x[2])


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
