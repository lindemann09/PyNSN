#!/usr/bin/env python
from expyriment import control, design, io, stimuli
from turning_knob import numerosity_production
import nosynum
from nosynum import expyriment_stimulus, pil_image

#control.set_develop_mode(True)

exp = design.Experiment(background_colour=(255,255,255))

control.initialize(exp)

maxnumber = 200

dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 400,
                       dot_diameter_mean=10,
                       dot_diameter_range=(2, 30),
                       dot_diameter_std=2,
                       dot_colour=(255, 57, 57))

make_process = None

control.start(exp, skip_ready_screen=True)
#######################
mouse = io.Mouse(show_cursor=False, track_button_events=True)


stimuli.TextLine("Creating stimuli....please wait").present()



stimuli.BlankScreen().present()

while(True):
    if make_process is None: # just the first time
        stimuli.TextLine("Please waiting...").present()
        da_sequence = None
    else:
        da_sequence = make_process.da_sequence
    max_da = nosynum.DotArray(n_dots=maxnumber, dot_array_definition=dot_array_def)
    make_process = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                                      method=nosynum.M_NO_FITTING,
                                                      save_images=True)
    if da_sequence is not None:
        #array_id_estim = design.randomize.rand_int(0, len(arrays) - 1)
        estim, rt, stim = numerosity_production(exp, mouse,
                                                da_sequence,
                                                maxnumber)
        print(estim, rt)

control.end()
