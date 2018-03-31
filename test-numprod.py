#!/usr/bin/env python
import time
from expyriment import control, design, io, stimuli, misc
import nosynum
from nosynum import expyriment_stimulus, pil_image
from turning_knob import numerosity_production

control.set_develop_mode(True)

cl = misc.Clock()

generator =  nosynum.RandomDotArrayGenerator(
                       stimulus_area_radius= 400,
                       dot_diameter_mean=10,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=1,
                       dot_colour=(230, 230, 230))


def prepare_processes(trials, forerun):
    processes = nosynum.ProcessContainer(forerun=forerun)
    for tr in trials:
        maxnumber = tr.get_factor("maxnumber")
        method = tr.get_factor("method")
        max_da = generator(n_dots=maxnumber)
        make_process = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                                          min_numerosity=20,
                                                          method=method,
                                                          area_colour=(0, 0, 0),
                                                          antialiasing=None,
                                                          save_images=False,
                                                          sqeeze_factor=.70)
        processes.add_process(make_process)

    while processes.data_avaiable < forerun:
        time.sleep(1)

    return processes


#########################################
exp = control.initialize()

bl = design.Block()
for m in [nosynum.M_TOTAL_CIRCUMFERENCE, nosynum.M_DENSITY, nosynum.M_NO_FITTING,\
            nosynum.M_CONVEX_HULL, nosynum.M_TOTAL_AREA]:
    tr = design.Trial()
    tr.set_factor("maxnumber", 200)
    tr.set_factor("method", m)
    for num in [20, 40, 60, 80, 100]:
        tr.set_factor("target", num)
        bl.add_trial(tr, copies=20)

bl.shuffle_trials()
exp.add_block(bl)


control.start(exp, skip_ready_screen=True)
#######################

stimuli.TextLine("Creating stimuli....please wait").present()

blank = stimuli.BlankScreen()
fixcross = stimuli.FixCross()
da_raw_file = io.OutputFile(suffix=".da-raw.csv",directory="data")
da_property_file = io.OutputFile(suffix=".da-prop.csv",directory="data")

exp.mouse.hide_cursor(track_button_events=True, track_motion_events=False)


processes = prepare_processes(exp.blocks[0].trials, forerun=1)
for tr_cnt, tr in enumerate(exp.blocks[0].trials):

    make_process = processes.pop_processes()
    if not make_process.data_available.is_set():
        stimuli.TextLine("Please wait...").present()

    da_sequence = make_process.da_sequence
    print(da_sequence.min_max_numerosity[1])

    if da_sequence is not None: # error handling
        blank.present()

        da_expy = expyriment_stimulus.ExpyrimentDASequence(da_sequence)
        da_expy.create_stimuli_from_pil_images()
        da_expy.preload(percent=60)

        fixcross.present()
        da_expy.preload(time=500, do_not_return_earlier=True)
        blank.present()
        da_expy.preload(time=500, do_not_return_earlier=True)
        stimuli.TextLine(str(234)).present()
        da_expy.preload(time=1000, do_not_return_earlier=True)
        blank.present()
        # cl.reset_stopwatch()
        da_expy.preload(time=500, do_not_return_earlier=True)
        da_expy.preload()  # ensure the rest is preloaded (should no be required)

        # print(cl.stopwatch_time)
        estim, rt, _, adjustments = numerosity_production(exp=exp, mouse=exp.mouse,
                                                expy_da_sequence=da_expy,
                                                start_value=da_sequence.min_max_numerosity[1])
        da_expy.unload()

        # saving data
        da_raw_file.write_line(da_sequence.get_csv(num_format="%6.0f", hash_column=True,
                                                   variable_names=(tr_cnt==0)))
        da_property_file.write_line(da_sequence.get_property_string(variable_names=(tr_cnt == 0)))

        # print(da_sequence.get_csv(num_format="%6.0f"))
        # print(da_expy.da_sequence.numerosity_correlations)

        exp.data.add([tr_cnt, tr.get_factor("target"), da_sequence.md5hash,
                      tr.get_factor("method"), estim, rt])
        if (estim is None):
            break
    else:
        print(da_sequence.error)


control.end()
