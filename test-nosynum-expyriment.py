#!/usr/bin/env python

import nosynum
from nosynum import expyriment_stimulus, pil_image

from expyriment import control, misc, stimuli

c = misc.Clock()
control.set_develop_mode(True)

exp = control.initialize()


dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=15,
                       dot_diameter_range=(2, 20),
                       dot_diameter_std=3)

max_da = nosynum.DotArray(n_dots=100, dot_array_definition=dot_array_def)
mp = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                        save_images=True,
                                        method=nosynum.M_DENSITY)
mp.join()

control.start()

blank = stimuli.BlankScreen()

dot_picture=None#"picts/pict1.png"

while(True):
    blank.present()
    # get next sequence and restart
    c.reset_stopwatch()
    print("DAS seq")
    expy_seq = expyriment_stimulus.ExpyrimentDASequence(da_sequence=mp.da_sequence)
    max_da = nosynum.DotArray(n_dots=100, dot_array_definition= dot_array_def)
    mp = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                       method=nosynum.M_DENSITY,
                                       save_images=True,
                                       auto_start_process=True)
    print(c.stopwatch_time)

    c.reset_stopwatch()
    print("EXPY")
    expy_seq.load_associated_images()
    #expy_seq.create_stimuli()
    #expy_seq.preload()
    print(c.stopwatch_time)

    expy_seq.stimuli[10].present()

    key, rt = exp.keyboard.wait()
    if key == misc.constants.K_RETURN:
        break

    expy_seq.unload()

#d.save_incremental_images(name="picts/xx", file_type="png", antialiasing=False, area_colour=(255,255,255))

control.end()
