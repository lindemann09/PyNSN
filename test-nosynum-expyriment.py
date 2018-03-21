#!/usr/bin/env python

import nosynum
from nosynum import expyriment_stimulus

from expyriment import control, misc

c = misc.Clock()
control.set_develop_mode(True)

exp = control.initialize()
control.start()

dot_diameter_size = 40
dot_diameter_range= (2, 50)
dot_diameter_std = 4
dot_picture=None#"picts/pict1.png"

while(True):
    c.reset_stopwatch()
    d = nosynum.DotArray(stimulus_area_radius= 300, n_dots=20,
            dot_diameter_range = dot_diameter_range,
            dot_diameter_mean = dot_diameter_size,
            dot_diameter_std = dot_diameter_std,
            dot_picture= dot_picture,
            dot_colour=(0,0,0))
    print(c.stopwatch_time)
    print "area",d.convex_hull_area
    c.reset_stopwatch()
    #print d.dots
    stim = expyriment_stimulus.create(d, area_colour=misc.constants.C_GREY, anti_aliasing=10)
    stim.preload()
    print(c.stopwatch_time)
    stim.present()

    key, rt = exp.keyboard.wait()
    if key == misc.constants.K_RETURN:
        break

#d.save_incremental_images(name="picts/xx", file_type="png", antialiasing=False, area_colour=(255,255,255))

control.end()
