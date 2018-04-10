#!/usr/bin/env python

import pynsn
from pynsn import expyriment_stimulus

from expyriment import control, misc, stimuli

control.defaults.window_size=(1100, 600)
c = misc.Clock()
control.set_develop_mode(True)
control.defaults.open_gl = False

exp = control.initialize()

squeeze = .7
generator = pynsn.DotArrayGenerator(
                       max_array_radius=300 * squeeze,
                       dot_diameter_mean=10,
                       dot_diameter_range=(5, 20),
                       dot_diameter_std=2,
                       minimum_gap=3)

def compare_stimulus(n_left, n_right,
                     pos_left=(-300,0),
                     pos_right=(300,0),
                     match_method=pynsn.DASequenceGenerator.NO_FITTING,
                     match_the_left=True):

    da_left = generator.make(n_dots=n_left)
    da_right= generator.make(n_dots=n_right)

    if match_the_left:
        b = da_right
        a = da_left
    else:
        a = da_right
        b = da_left

    if match_method == pynsn.DASequenceGenerator.DENSITY50_50:
        a._match_density(density=b.density, adaptation_CH2TA_ratio=0.5)
    elif match_method == pynsn.DASequenceGenerator.CONVEX_HULL:
        a._match_convex_hull_area(convex_hull_area=b.convex_hull_area)
    elif match_method == pynsn.DASequenceGenerator.MEAN_DIAMETER:
        a._match_mean_dot_diameter(b.mean_dot_diameter)
    elif match_method == pynsn.DASequenceGenerator.TOTAL_AREA:
        a.match_total_area(total_area=b.total_area)
    elif match_method == pynsn.DASequenceGenerator.TOTAL_CIRCUMFERENCE:
        a._match_total_circumference(total_circumference=b.total_circumference)
    elif match_method == pynsn.DASequenceGenerator.NO_FITTING:
        pass
    else:
        raise Warning("Unknown method {}. Using NO_FITTING.".format(match_method))


    ok, mesg = a.realign()

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
                        match_method=pynsn.DASequenceGenerator.NO_FITTING,
                        match_the_left=True)

stim.present()

key, rt = exp.keyboard.wait()
control.end()
