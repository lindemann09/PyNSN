#!/usr/bin/env python
import numpy as np
from time import sleep
from expyriment import control, design, io, stimuli, misc
from turning_knob import numerosity_production
import nosynum
from nosynum import expyriment_stimulus, pil_image

control.set_develop_mode(True)

exp = design.Experiment()#background_colour=(255,255,255))

control.initialize(exp)

maxnumber = 200

dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=10,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=1,
                       dot_colour=(230, 230, 230))

make_process = None

control.start(exp, skip_ready_screen=True)
#######################
mouse = io.Mouse(show_cursor=False, track_button_events=True, track_motion_events=False)


stimuli.TextLine("Creating stimuli....please wait").present()

cl = misc.Clock()
blank = stimuli.BlankScreen()
fixcross = stimuli.FixCross()

da_sequence = None

forerun = 1

tmp = nosynum.ProcessContainer(forerun=forerun)
for x in range(20):
    max_da = nosynum.DotArray(n_dots=maxnumber, dot_array_definition=dot_array_def)
    make_process = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                                      method=nosynum.M_TOTAL_CIRCUMFERENCE,
                                                      area_colour=(0,0,0),
                                                      antialiasing=None,
                                                      save_images=False)
    tmp.add_process(make_process)

while tmp.data_avaiable<forerun:
    sleep(1)
processes_A = tmp


while(True):

    make_process = processes_A.pop_processes()

    if not make_process.data_available.is_set():
        stimuli.TextLine("Please wait...").present()

    da_sequence = make_process.da_sequence
    print(da_sequence.error)
    if da_sequence is not None:
        blank.present()

        da_expy = expyriment_stimulus.ExpyrimentDASequence(da_sequence)
        da_expy.create_stimuli_from_pil_images()

        da_expy.preload(time=500)
        fixcross.present()
        da_expy.preload(time=200)
        blank.present()
        da_expy.preload(time=500)
        stimuli.TextLine(str(234)).present()
        da_expy.preload(time=1000, do_not_return_earlier=True)
        blank.present()
        da_expy.preload(time=300, do_not_return_earlier=True)
        da_expy.preload() # the rest

        print(da_expy.da_sequence.get_csv(num_format="%6.0f"))
        print(da_expy.da_sequence.property_string)
        print(da_expy.da_sequence.numerosity_correlations)
        print(da_expy.da_sequence.variances)
        estim, rt, stim = numerosity_production(exp=exp, mouse=mouse, expy_da_sequence=da_expy,
                                                start_value=da_sequence.min_max_numerosity[1])
        blank.present()
        da_expy.unload()
        print(estim, rt)
        if (estim is None):
            break

control.end()
