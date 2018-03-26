#!/usr/bin/env python

from PIL import Image
import nosynum
from nosynum import expyriment_stimulus, pil_image

from expyriment import control, misc, stimuli
import pygame

c = misc.Clock()
control.set_develop_mode(True)

exp = control.initialize()


dot_array_def = nosynum.DotArrayDefinition(
                       stimulus_area_radius= 300,
                       dot_diameter_mean=5,
                       dot_diameter_range=(2, 10),
                       dot_diameter_std=2)

max_da = nosynum.DotArray(n_dots=100, dot_array_definition=dot_array_def)
mp = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                        save_images=False,
                                        method=nosynum.M_NO_FITTING)
mp.join()

control.start(skip_ready_screen=True)

blank = stimuli.BlankScreen()

dot_picture=None#"picts/pict1.png"

while(True):
    blank.present()
    # get next sequence and restart
    c.reset_stopwatch()
    print("DAS seq")
    da_sequence = mp.da_sequence
    expy_seq = expyriment_stimulus.ExpyrimentDASequence(da_sequence=da_sequence)
    max_da = nosynum.DotArray(n_dots=100, dot_array_definition= dot_array_def)
    mp = pil_image.PILMakeDASequenceProcess(max_dot_array=max_da,
                                       method=nosynum.M_DENSITY,
                                       save_images=False)
    mp.start()

    print(c.stopwatch_time)
    c.reset_stopwatch()
    print("EXPY")
    expy_seq.load_associated_images()
    expy_seq.preload()

    print(c.stopwatch_time)
    c.reset_stopwatch()

    mp.start()
    expy_seq.stimuli[10].present()

    key, rt = exp.keyboard.wait()
    if key == misc.constants.K_RETURN:
        break

    expy_seq.unload()

#d.save_incremental_images(name="picts/xx", file_type="png", antialiasing=False, area_colour=(255,255,255))

control.end()
