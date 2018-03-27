'''
Created on Nov 4, 2010

@author: oliver
'''
import sys
import time

from expyriment import stimuli
from expyriment.misc import constants, Clock


class Knob:
    large_step = 1
    def __init__(self, min_value, max_value, time_interval, start_value=0):
        self.position = None
        self.value = None
        self._min_value = min_value
        self._max_value = max_value
        self._time_interval = time_interval
        self._last_change = 0
        if sys.platform == 'win32':
            self._cpu_time = time.clock
        else:
            self._cpu_time = time.time
        self.reset(start_value=start_value)

    def reset(self, start_value=0):
        self.position = 0
        self.value = start_value

    def turn_clockwise(self):
        self.position += 1
        self._change_value(True)

    def turn_counterclockwise(self):
        self.position -= 1
        self._change_value(False)

    def _change_value(self, increase_value=True):
        if (self._cpu_time()*1000) - self._last_change >= self._time_interval:
            step = 1
        else:
            step = self.large_step
        self._last_change = self._cpu_time()*1000
        if increase_value:
            self.value += step
        else:
            self.value -= step

        if self.value < self._min_value:
            self.value = self._min_value
        elif self.value > self._max_value:
            self.value = self._max_value


def  numerosity_production(exp, mouse, expy_da_sequence, preload=False, background=None,
                          start_value=0):
    """numerosity production

    expy_da_sequence: associated stimuli have to be created

    returns knob value, rt, and last stimuli, adjustments (all actions (value, time))
    """
    if background is None:
        background = stimuli.BlankScreen()

    background.present()
    if  preload:
        expy_da_sequence.preload()

    knob = Knob(expy_da_sequence.da_sequence.min_max_numerosity[0],
                expy_da_sequence.da_sequence.min_max_numerosity[1],
                time_interval=0, start_value=start_value)

    exp.clock.reset_stopwatch()
    cl = Clock()
    adjustments =[]

    goOn = True
    changed = True
    stim = None
    while goOn:
        if changed:
            stim = expy_da_sequence.get_stimulus_numerosity(knob.value)
            stim.present(clear=False, update=True)
            adjustments.append((knob.value, exp.clock.stopwatch_time))
            changed = False

        k = exp.keyboard.check()
        m = mouse.get_last_button_down_event()
        if k == constants.K_LEFT or m == 3:
            knob.turn_counterclockwise()
            #tone2.play()
            changed = True
        elif k == constants.K_RIGHT or m == 4:
            knob.turn_clockwise()
            #tone.play()
            changed = True
        elif (k == constants.K_SPACE or m == 0) and knob.value>0:
            goOn = False
        elif (k == constants.K_q):
            goOn = False
            knob.value = None
        else:
            time.sleep(0.05)
            if cl.stopwatch_time>3000:
                mouse.clear()
                cl.reset_stopwatch()

    rt = exp.clock.stopwatch_time
    background.present()
    if preload:
        expy_da_sequence.unload()
    return knob.value, rt, stim, adjustments
