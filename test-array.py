#!/usr/bin/env python
from __future__ import absolute_import, print_function, division

from expyriment import misc, control
import pynsn
from pynsn import expyriment_stimulus, DotArrayGenerator, GeneratorLogger

if __name__ == "__main__":
    cl = misc.Clock()

    logger = GeneratorLogger(log_filename="log/test", override_log_files=True,
                             log_colours=True, properties_different_colour=True)
    generator = DotArrayGenerator(
        target_area_radius=300,
        item_diameter_mean=20,
        item_diameter_range=(5, 40),
        item_diameter_std=2,
        item_colour="skyblue",
        minimum_gap=5)

    reference = generator.make(n_dots=25, logger=None)
    reference.center_array()
    reference.realign()
    print(reference.get_features_text(with_object_id=False, extended_format=False ))

    reference2 = reference.copy()
    reference2.match(match_features=pynsn.Coverage(0.10, match_ratio_fieldarea2totalarea=1), center_array=True)
    #reference2.realign()
    print(reference2.get_features_text(with_object_id=False, extended_format=False))

    # ds = make_dot_array_sequence(reference_dot_array=reference,
    #                              logger=logger,
    #                              extra_space=100,
    #                              match_properties=[pynsn.Coverage()],
    #                              # match_methods=[DASequenceGenerator.CONVEX_HULL,
    #                              #              DASequenceGenerator.DENSITY_ONLY_AREA],
    #                              min_max_numerosity=[10, 30])
    #
    # prop = ds.get_features_dict()
    # exit()
    # # print(ds.get_numerosity_correlations())
    # # print(ds.variances)
    # # print(ds.numerosity_correlations)
    #
    #
    #
    # if True:
    #
    #     a = generator.make(n_dots=20, logger=logger)
    #     b = generator.make(n_dots=40, logger=logger)
    #     b.features.change(colour="blue")
    #     # b.match(mean_dot_diameter=83, convex_hull_area=523)
    #     x = a.copy()
    #     x.join(b)
    #     logger.log(x)
    #
    #
    # else:
    #     p = DASequenceMakeProcess(target_dot_array=reference,
    #                               match_method=[DASequenceGenerator.CONVEX_HULL,
    #                                                  DASequenceGenerator.DENSITY_ONLY_AREA],
    #                               min_numerosity=10, logger=logger,
    #                               extra_space=100)
    #     p.start()
    #     p1 = DASequenceMakeProcess(target_dot_array=reference, match_method=DASequenceGenerator.CONVEX_HULL,
    #                                min_numerosity=10, logger=logger,
    #                                extra_space=100)
    #     p1.start()
    #     x = p.da_sequence

    # exit()

    # print(x.get_csv(colour_column=True))
    control.set_develop_mode(True)
    exp = control.initialize()
    control.start(exp, skip_ready_screen=True)

    stimx = expyriment_stimulus.ExprimentDotArray(reference,
                                                  colour_target_area="dimgray",
                                                  # colour_convex_hull_dots=(255, 0, 0),
                                                  # colour_convex_hull_positions=(255, 200, 0),
                                                  # colour_center_of_mass=(255, 0, 0),
                                                  #  colour_center_of_outer_positions=(0,0,200),
                                                  antialiasing=True)

    stimx.present()
    exp.keyboard.wait()
    stimx = expyriment_stimulus.ExprimentDotArray(reference2,
                                                  colour_target_area="dimgray",
                                                  # colour_convex_hull_dots=(255, 0, 0),
                                                  # colour_convex_hull_positions=(255, 200, 0),
                                                  # colour_center_of_mass=(255, 0, 0),
                                                  #  colour_center_of_outer_positions=(0,0,200),
                                                  antialiasing=True)
    stimx.present()
    exp.keyboard.wait()

    control.end()
