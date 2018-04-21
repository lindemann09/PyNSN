#!/usr/bin/env python
import time
from expyriment import control, design, io, stimuli, misc
import pynsn
from pynsn import continuous_property as cp
from pynsn import expyriment_stimulus, pil_image, DASequenceGenerator
from turning_knob import numerosity_production

control.set_develop_mode(True)

cl = misc.Clock()

generator =  pynsn.DotArrayGenerator(
                       max_array_radius= 200,
                       dot_diameter_mean=10,
                       dot_diameter_range=(5, 15),
                       dot_diameter_std=1,
                       dot_colour= None,
                       logger= pynsn.GeneratorLogger(log_filename="log/test", override_log_files=True)
                )

#########################################
exp = control.initialize()

bl = design.Block()
for m in [cp.Circumference(),
    DASequenceGenerator.TOTAL_CIRCUMFERENCE,
          DASequenceGenerator.DENSITY50_50,
          DASequenceGenerator.NO_FITTING,
          DASequenceGenerator.CONVEX_HULL,
          DASequenceGenerator.TOTAL_AREA]:
    tr = design.Trial()
    tr.set_factor("maxnumber", 100)
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

exp.mouse.hide_cursor(track_button_events=True, track_motion_events=False)

def prepare_processes(trials):
    processes = []
    for tr in trials:
        maxnumber = tr.get_factor("maxnumber")
        method = tr.get_factor("method")
        max_da = generator.make(n_dots=maxnumber, inhibit_logging=True)
        make_process = pynsn.DASequenceGeneratorProcess(max_dot_array=max_da,
                                                          min_numerosity=20,
                                                          match_methods=method,
                                                          extra_space=100,
                                                          logger=generator.logger)
        processes.append(make_process)
    return processes

def make_da_sequence(trial):
    maxnumber = trial.get_factor("maxnumber")
    method = trial.get_factor("method")
    max_da = generator.make(n_dots=maxnumber, inhibit_logging=True)
    da_gen = pynsn.DASequenceGenerator(max_dot_array=max_da,
                                        logger=generator.logger)
    return da_gen.make(min_numerosity=20, match_properties=method, extra_space=100)


processes = prepare_processes(bl.trials)

for tr_cnt, tr in enumerate(exp.blocks[0].trials):

    if True: # process
        make_process = processes.pop()
        make_process.start()
        if not make_process.data_available.is_set():
            stimuli.TextLine("Please wait...").present()
        da_sequence = make_process.da_sequence
    else:
        da_sequence = make_da_sequence(tr)

    print(da_sequence.min_max_numerosity[1])

    if da_sequence is not None: # error handling
        blank.present()

        da_expy = expyriment_stimulus.ExpyrimentDASequence(da_sequence)
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
        # print(da_sequence.get_csv(num_format="%6.0f"))
        # print(da_expy.da_sequence.numerosity_correlations)

        exp.data.add([tr_cnt, tr.get_factor("target"), da_sequence.object_id,
                      tr.get_factor("method"), estim, rt])
        if (estim is None):
            break
    else:
        print(da_sequence.error)


control.end()
