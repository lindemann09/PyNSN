#!/usr/bin/env python

from PIL import Image
import nosynum
from nosynum import expyriment_stimulus, pil_image

from expyriment import control, misc, stimuli
import pygame

c = misc.Clock()
control.set_develop_mode(True)

exp = control.initialize()


squeeze = .7
dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 300*squeeze,
                       dot_diameter_mean=5,
                       dot_diameter_range=(2, 10),
                       dot_diameter_std=2)

def compare_stimulus(n_left, n_right, pos_left=(-300,0), pos_right=(300,0),
                     match_method=nosynum.M_NO_FITTING,
                     match_the_left=True):

    da_left = nosynum.DotArray(n_dots=n_left, dot_array_definition=dot_array_def)

    da_right= nosynum.DotArray(n_dots=n_right, dot_array_definition=dot_array_def)

    if match_the_left:
        b = da_right
        a = da_left
    else:
        a = da_right
        b = da_left

    if match_method == nosynum.M_DENSITY:
        a.fit_density(density=b.density, ratio_area_convex_hull_adaptation=0.5)
    elif match_method == nosynum.M_CONVEX_HULL:
        a.fit_convex_hull_area(convex_hull_area=b.convex_hull_area)
    elif match_method == nosynum.M_ITEM_SIZE:
        a.fit_mean_item_size(b.mean_dot_diameter)
    elif match_method == nosynum.M_TOTAL_AREA:
        a.fit_total_area(total_area=b.total_area)
    elif match_method == nosynum.M_TOTAL_CIRCUMFERENCE:
        a.fit_total_circumference(total_circumference=b.total_circumference)
    elif match_method == nosynum.M_NO_FITTING:
        pass
    else:
        raise Warning("Unknown method {}. Using NO_FITTING.".format(match_method))

    # stimuli
    left = expyriment_stimulus.ExprimentPILImage(pil_image.create(dot_array=da_left),
                                              position=pos_left)
    right = expyriment_stimulus.ExprimentPILImage(pil_image.create(dot_array=da_right),
                                              position=pos_right)
    stim = stimuli.BlankScreen()
    left.plot(stim)
    right.plot(stim)
    return stim

control.start(skip_ready_screen=True)

c.reset_stopwatch()

stim = compare_stimulus(50, 120)

print(c.stopwatch_time)
c.reset_stopwatch()
stim.preload()
print(c.stopwatch_time)

stim.present()

key, rt = exp.keyboard.wait()
control.end()
