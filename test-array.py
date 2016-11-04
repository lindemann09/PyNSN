#!/usr/bin/env python

import math
import random
import numpy as np
from nosynum import *
from expyriment import misc
c = misc.Clock()

dot_diameter_range=(2, 20)
number_list = range(10, 200)

print "creating "
dal = DotArrayList(number_list, subfolder="picts")
dal.create_method(method=3, dot_diameter = 20, dot_diameter_std=2)

#dal2 = DotArrayList(number_list, subfolder="picts")
#dal2.create_method(method=2, dot_diameter = 10)
#dal.extend(dal2)
#dal2.create_method(method=3, dot_diameter = 10)
#dal.extend(dal2)

dal.analysis()
dal.save_arrays_csv_text(filename="dot_arrays.csv")
dal.save_properties("array_info.csv")

#print "saving"
#dal.save_images(name="test", file_type="JPEG", antialiasing=False)

#import matplotlib.pyplot as plt
#plt.plot(data[:,0], data[:, 4])
#plt.show()
