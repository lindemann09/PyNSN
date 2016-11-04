#!/usr/bin/env python

from nosynum import DotArray
from expyriment import control, misc

control.set_develop_mode(True)



exp = control.initialize()
control.start()

dot_diameter_size = None

dot_diameter_range=(2, 20)
dot_diameter_std = 4

while(True):
    if dot_diameter_size == 8:
        dot_diameter_size = 15
    else:
        dot_diameter_size = 8

    d = DotArray(stimulus_area_radius= 250, n_dots=50,
            dot_diameter_range = dot_diameter_range,
            dot_diameter_mean = dot_diameter_size,
            dot_diameter_std = dot_diameter_std)
    print "area",d.convex_hull_area
    #print d.dots
    d.create_expyriment_stimulus(area_colour=(50,50,50)).present()
    key, rt = exp.keyboard.wait()
    if key == misc.constants.K_RETURN:
        break

control.end()
