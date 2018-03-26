#!/usr/bin/env python
import numpy as np
from time import sleep
from expyriment import control, design, io, stimuli, misc
from turning_knob import numerosity_production
import nosynum
from nosynum import expyriment_stimulus, pil_image
control.set_develop_mode(True)


def numprod(da_sequence):

    blank = stimuli.BlankScreen()
    fixcross = stimuli.FixCross()

    if da_sequence is not None:
        blank.present()

        da_expy = expyriment_stimulus.ExpyrimentDASequence(da_sequence)
        da_expy.create_stimuli_from_pil_images()

        da_expy.preload(percent=60)
        fixcross.present()
        da_expy.preload(time=500)
        blank.present()
        da_expy.preload(time=500)
        stimuli.TextLine(str(234)).present()
        da_expy.preload(time=1000, do_not_return_earlier=True)
        blank.present()
        # cl.reset_stopwatch()
        da_expy.preload(time=500, do_not_return_earlier=True)
        da_expy.preload()  # the rest

        # print(cl.stopwatch_time)
        estim, rt, stim = numerosity_production(exp=exp, mouse=mouse, expy_da_sequence=da_expy,
                                                start_value=da_sequence.min_max_numerosity[1])
        blank.present()
        da_expy.unload()
        return estim, rt


exp = design.Experiment()#background_colour=(255,255,255))

control.initialize(exp)

maxnumber = 200


dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 400,
                       dot_diameter_mean=10,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=1,
                       dot_colour=(230, 230, 230))

make_process = None

control.start(exp, skip_ready_screen=True)
#######################
mouse = io.Mouse(show_cursor=False, track_button_events=True, track_motion_events=False)

da_raw_file = io.OutputFile(suffix=".da-raw.csv",directory="data")
da_property_file = io.OutputFile(suffix=".da-prop.csv",directory="data")

stimuli.TextLine("Creating stimuli....please wait").present()

cl = misc.Clock()


methods = [nosynum.M_TOTAL_CIRCUMFERENCE, nosynum.M_DENSITY, nosynum.M_NO_FITTING,\
            nosynum.M_CONVEX_HULL, nosynum.M_TOTAL_AREA]*20
design.randomize.shuffle_list(methods)

print(methods)
da_sequence = None
forerun = 3
processes = nosynum.ProcessContainer(forerun=forerun)
for m in methods:
    max_da = nosynum.DotArray(n_dots=maxnumber, dot_array_definition=dot_array_def)
    make_process = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                                      method=m,
                                                      area_colour=(0,0,0),
                                                      antialiasing=None,
                                                      save_images=False,
                                                      sqeeze_factor=.70)
    processes.add_process(make_process)


while processes.data_avaiable<forerun:
    sleep(1)

write_varnames = True
while(True):
    make_process = processes.pop_processes()

    if not make_process.data_available.is_set():
        stimuli.TextLine("Please wait...").present()

    da_sequence = make_process.da_sequence
    if da_sequence.error is not None:
        print(da_sequence.error)
    # print(da_sequence.method)

    if da_sequence is not None:
        estim, rt = numprod(da_sequence)

        da_raw_file.write_line(da_sequence.get_csv(num_format="%6.0f", hash_column=True,
                                               variable_names=write_varnames))
        da_property_file.write_line(da_sequence.get_property_string(varnames=True))
        write_varnames = False
        # print(da_sequence.get_csv(num_format="%6.0f"))
        # print(da_expy.da_sequence.numerosity_correlations)

        if (estim is None):
            break

control.end()
