#!/usr/bin/env python
import numpy as np
from expyriment import control, design, io, stimuli, misc
from turning_knob import numerosity_production
import nosynum
from nosynum import expyriment_stimulus, pil_image

control.set_develop_mode(True)

exp = design.Experiment()#background_colour=(255,255,255))

control.initialize(exp)

maxnumber = 300

dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 200,
                       dot_diameter_mean=10,
                       dot_diameter_range=(2, 30),
                       dot_diameter_std=2,
                       dot_colour=(255, 57, 57))

make_process = None

control.start(exp, skip_ready_screen=True)
#######################
mouse = io.Mouse(show_cursor=False, track_button_events=True, track_motion_events=False)


stimuli.TextLine("Creating stimuli....please wait").present()

cl = misc.Clock()
blank = stimuli.BlankScreen()
fixcross = stimuli.FixCross()

da_sequence = None

while(True):

    if make_process is not None:
        if not make_process.sequence_available.is_set():
            stimuli.TextLine("Please wait...").present()
            make_process.join()
        da_sequence = make_process.da_sequence

    max_da = nosynum.DotArray(n_dots=maxnumber, dot_array_definition=dot_array_def)
    make_process = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                                      method=nosynum.M_NO_FITTING,
                                                      area_colour=(0,0,0),
                                                      save_images=True,
                                                      auto_start_process=False)

    if da_sequence is not None:
        blank.present()
        da_expy = expyriment_stimulus.ExpyrimentDASequence(da_sequence)
        da_expy.load_associated_images()
        da_expy.preload(time=1000)

        fixcross.present()
        da_expy.preload(time=500)
        blank.present()
        da_expy.preload(time=500)
        stimuli.TextLine(str(234)).present()
        da_expy.preload(time=1000, do_not_return_earlier=True)
        blank.present()
        da_expy.preload(time=500, do_not_return_earlier=True)
        da_expy.preload() # the rest
        make_process.start()
        estim, rt, stim = numerosity_production(exp=exp, mouse=mouse, expy_da_sequence=da_expy,
                                                start_value=da_sequence.min_max_numerosity[1])
        blank.present()
        da_expy.unload()
        print(estim, rt)
        if (estim is None):
            break
    else:
        make_process.start()

control.end()
