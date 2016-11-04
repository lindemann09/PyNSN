#!/usr/bin/env python

from nosynum import DotArray
from expyriment import control, misc

control.set_develop_mode(True)

exp = control.initialize()
control.start()

dot_diameter_size = 110
dot_diameter_range= None #(2, 20)
dot_diameter_std = 4
dot_picture=None#"picts/pict1.png"

while(True):
    d = DotArray(stimulus_area_radius= 400, n_dots=20,
            dot_diameter_range = dot_diameter_range,
            dot_diameter_mean = dot_diameter_size,
            dot_diameter_std = dot_diameter_std,
            dot_picture= dot_picture,
            dot_colour=(0,0,0))
    print "area",d.convex_hull_area
    #print d.dots
    d.create_expyriment_stimulus(area_colour=(255,255,255)).present()
    key, rt = exp.keyboard.wait()
    if key == misc.constants.K_RETURN:
        break

d.save_incremental_images(name="picts/xx",
                         file_type="png",
                         antialiasing=False,
                          area_colour=(255,255,255))

control.end()
