#!/usr/bin/env python
from __future__ import absolute_import, print_function, division

from expyriment import misc, control
import pynsn
from pynsn import expyriment_stimulus, DotArrayGenerator, \
    DASequenceGeneratorProcess, DASequenceGenerator, GeneratorLogger, DotArrayProperties



if __name__ == "__main__":
    cl = misc.Clock()

    maxnumber = 20

    logger = GeneratorLogger(log_filename="log/test", override_log_files=True,
                             log_colours=True, properties_different_colour=True)
    generator = DotArrayGenerator(
        max_array_radius=300,
        dot_diameter_mean=25,
        dot_diameter_range=(5, 40),
        dot_diameter_std=2,
        dot_colour="skyblue",
        minimum_gap=5,
        logger=logger)
    max_da = generator.make(n_dots=maxnumber, inhibit_logging=False)

    g = DASequenceGenerator(max_dot_array=max_da,  logger=logger)
    ds = g.make(extra_space=100, match_properties=None,
                #match_methods=[DASequenceGenerator.CONVEX_HULL,
                #              DASequenceGenerator.DENSITY_ONLY_AREA],
                min_numerosity=10)

    prop = ds.get_properties()
    #print(ds.get_numerosity_correlations())
    #print(ds.variances)
    #print(ds.numerosity_correlations)



    if True:

        a = generator.make(n_dots=20, inhibit_logging=True)
        b = generator.make(n_dots=40, inhibit_logging=True)
        b.features.change(colour="blue")
        # b.match(mean_dot_diameter=83, convex_hull_area=523)
        x = a.copy()
        x.join(b)
        logger.log(x)
        print(x.get_properties_split_by_colours().get_csv())


    else:
        p = DASequenceGeneratorProcess(max_dot_array=max_da,
                                       match_method=[DASequenceGenerator.CONVEX_HULL,
                                                     DASequenceGenerator.DENSITY_ONLY_AREA],
                                       min_numerosity=10, logger=logger,
                                       extra_space=100)
        p.start()
        p1 = DASequenceGeneratorProcess(max_dot_array=max_da, match_method=DASequenceGenerator.CONVEX_HULL,
                                        min_numerosity=10, logger=logger,
                                        extra_space=100)
        p1.start()
        x = p.da_sequence

    # exit()

    # print(x.get_csv(colour_column=True))
    control.set_develop_mode(True)
    exp = control.initialize()
    control.start(exp, skip_ready_screen=True)

    stimx = expyriment_stimulus.ExprimentDotArray(x,
                                                 colour_area="dimgray",
                                                 # colour_convex_hull_dots=(255, 0, 0),
                                                 # colour_convex_hull_positions=(255, 200, 0),
                                                 # colour_center_of_mass=(255, 0, 0),
                                                 #  colour_center_of_outer_positions=(0,0,200),
                                                 antialiasing=True)

    stimx.present()
    exp.keyboard.wait()


    stima = expyriment_stimulus.ExprimentDotArray(a, colour_area="dimgray", antialiasing=True)
    stimb = expyriment_stimulus.ExprimentDotArray(b, colour_area="dimgray",antialiasing=True)
    stimx.preload()
    stima.preload()
    stimb.preload()
    while True:
        stima.present()
        exp.keyboard.wait()
        stimx.present()
        exp.keyboard.wait()
        stimb.present()
        exp.keyboard.wait()
        stimx.present()
        exp.keyboard.wait()

    control.end()
