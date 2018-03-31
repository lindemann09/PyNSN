#!/usr/bin/env python

import nosynum
from nosynum import expyriment_stimulus, pil_image

from expyriment import control, misc, stimuli

control.defaults.window_size=(1100, 600)
c = misc.Clock()
control.set_develop_mode(True)
control.defaults.open_gl = False

exp = control.initialize()

squeeze = .7
generator = nosynum.RandomDotArrayGenerator(
                       stimulus_area_radius= 300*squeeze,
                       dot_diameter_mean=10,
                       dot_diameter_range=(5, 20),
                       dot_diameter_std=2,
                       minimum_gap=3)

def compare_stimulus(n_left, n_right,
                     pos_left=(-300,0),
                     pos_right=(300,0),
                     match_method=nosynum.M_NO_FITTING,
                     match_the_left=True):

    da_left = generator.make(n_dots=n_left)
    da_right= generator.make(n_dots=n_right)

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
    elif match_method == nosynum.M_MEAN_DIAMETER:
        a.fit_mean_dot_diameter(b.mean_dot_diameter)
    elif match_method == nosynum.M_TOTAL_AREA:
        a.fit_total_area(total_area=b.total_area)
    elif match_method == nosynum.M_TOTAL_CIRCUMFERENCE:
        a.fit_total_circumference(total_circumference=b.total_circumference)
    elif match_method == nosynum.M_NO_FITTING:
        pass
    else:
        raise Warning("Unknown method {}. Using NO_FITTING.".format(match_method))


    ok = a.realign(minimum_gap=generator.minimum_gap) # todo: check of error

    # stimuli
    left = expyriment_stimulus.ExprimentDotArray( dot_array=da_left,  position=pos_left)
    right = expyriment_stimulus.ExprimentDotArray(dot_array=da_right, position=pos_right)

    stim = stimuli.BlankScreen()
    left.plot(stim)
    right.plot(stim)
    return stim




control.start(skip_ready_screen=True)

c.reset_stopwatch()

stim = compare_stimulus(20, 120,
                        match_method=nosynum.M_NO_FITTING,
                        match_the_left=True)

stim.present()

key, rt = exp.keyboard.wait()
control.end()
