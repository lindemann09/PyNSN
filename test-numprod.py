#!/usr/bin/env python
from expyriment import control, design, io, stimuli
from turning_knob import numerosity_production
from nosynum import get_list_of_incremental_dot_arrays

#control.set_develop_mode(True)

exp = design.Experiment(background_colour=(255,255,255))

control.initialize(exp)

control.start(exp, skip_ready_screen=True)
#######################
mouse = io.Mouse(show_cursor=False, track_button_events=True)

maxnumber = 20

stimuli.TextLine("Creating stimuli....please wait").present()
arrays = get_list_of_incremental_dot_arrays(max_n_dots = maxnumber,
                                            area_colour = (255, 255, 255),
                                            area_radius=400,
                                            dot_size = 110,
                                            dot_picture="picts/pict2.png",
                                            position=(0, 0),
                                            background_colour_pil="white")

stimuli.BlankScreen().present()

while(True):
    array_id_estim = design.randomize.rand_int(0, len(arrays) - 1)
    estim, rt, stim = numerosity_production(exp, mouse,
                        arrays[array_id_estim], maxnumber)
    print(estim, rt)

control.end()
