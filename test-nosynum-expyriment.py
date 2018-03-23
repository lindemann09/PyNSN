#!/usr/bin/env python

import nosynum
from nosynum import expyriment_stimulus, pil_image

from expyriment import control, misc, stimuli

c = misc.Clock()
control.set_develop_mode(True)

exp = control.initialize()
control.start()

blank = stimuli.BlankScreen()

dot_picture=None#"picts/pict1.png"

dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=15,
                       dot_diameter_range=(2, 20),
                       dot_diameter_std=3)

max_da = nosynum.DotArray(n_dots=100, dot_array_definition=dot_array_def)
mp = nosynum.MakeDASequenceProcess(max_dot_array=max_da,
                                   method=nosynum.M_DENSITY,
                                   auto_start_process=True)

while(True):
    blank.present()

    # get next sequence and restart
    seq = mp.da_sequence
    max_da = nosynum.DotArray(n_dots=100, dot_array_definition= dot_array_def)
    mp = nosynum.MakeDASequenceProcess(max_dot_array=max_da,
                                       method=nosynum.M_DENSITY,
                                       auto_start_process=False)

    mp.start()

    c.reset_stopwatch()
    print("pil")
    pil_image.dasequence2images(dot_array_sequence=seq, image_filename="asdf")
    print(c.stopwatch_time)
    c.reset_stopwatch()
    print("EXPY")
    expy_seq = expyriment_stimulus.ExpyrimentDASequence(da_sequence=seq)
    expy_seq.create_stimuli(area_colour=misc.constants.C_GREY, anti_aliasing=0)
    print(c.stopwatch_time)

    expy_seq.preload()
    expy_seq.stimuli[10].present()

    key, rt = exp.keyboard.wait()
    if key == misc.constants.K_RETURN:
        break

    expy_seq.unload()

#d.save_incremental_images(name="picts/xx", file_type="png", antialiasing=False, area_colour=(255,255,255))

control.end()
