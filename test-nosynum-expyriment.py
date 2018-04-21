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
                     match_methods=None,
                     match_the_left=True):

    da_left = generator.make(n_dots=n_left)
    da_right= generator.make(n_dots=n_right)
    da_right.features.change(colour="red")
    if match_the_left:
        b = da_right
        a = da_left
    else:
        a = da_right
        b = da_left

    if match_methods is not None:
        a.match(match_methods)

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
                        match_methods=[],
                        match_the_left=True)

stim.present()

key, rt = exp.keyboard.wait()
control.end()
