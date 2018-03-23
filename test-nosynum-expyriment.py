#!/usr/bin/env python

import nosynum
from nosynum import expyriment_stimulus

from expyriment import control, misc

c = misc.Clock()
control.set_develop_mode(True)

exp = control.initialize()
control.start()


dot_picture=None#"picts/pict1.png"

dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=15,
                       dot_diameter_range=(2, 20),
                       dot_diameter_std=3)

while(True):
    # get next sequence and restart
    c.reset_stopwatch()
    print("getting sequence and restart")
    try:
        seq = mp.da_sequence
    except:
        seq = None
    max_da = nosynum.DotArray(n_dots=100, dot_array_definition= dot_array_def)
    mp = nosynum.MakeDASequenceProcess(max_dot_array=max_da,
                                       method=nosynum.sequence_methods.DENSITY)
    mp.sequence_available.wait()
    print(c.stopwatch_time)
    c.reset_stopwatch()

    if seq is not None:
        print("making stimuli")
        expy_seq = expyriment_stimulus.ExpyrimentDASequence(da_sequence=mp.da_sequence)
        expy_seq.create_stimuli(area_colour=misc.constants.C_GREY, anti_aliasing=10)
        expy_seq.preload()

        print(c.stopwatch_time)
        c.reset_stopwatch()

        expy_seq.stimuli[10].present()

        key, rt = exp.keyboard.wait()
        if key == misc.constants.K_RETURN:
            break

        expy_seq.unload()
#d.save_incremental_images(name="picts/xx", file_type="png", antialiasing=False, area_colour=(255,255,255))

control.end()
