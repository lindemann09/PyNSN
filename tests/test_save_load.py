import os
import unittest
import tempfile

import pynsn
from pynsn import NSNFactory, ImageColours
from pynsn import distributions as distr
from pynsn import fit

TEMP_FLD = os.path.join(tempfile.gettempdir(), "pynsn_unit_test")
print("Unittest folder {}".format(TEMP_FLD))

#FIXME not pointarray tests so far


class TestArrays(unittest.TestCase):

    def setUp(self):
        N = 20

        self.my_colours = ImageColours(target_area="#EEEEEE",
                                       background=None,
                                       opacity_object=0.9,
                                       default_object_colour="darkmagenta",
                                       field_area_positions="magenta",
                                       field_area="gray",
                                       center_of_field_area="red",
                                       center_of_mass="magenta"
                                       )  # FIXME test opacity

        factory = NSNFactory(target_area_radius=200)
        factory.set_appearance_dot(diameter=(40, 10, 30),
                                   attributes=distr.Levels(["blue", "green"],
                                                           exact_weighting=True))
        self.dot_stim = factory.create_random_array(n_objects=N)

        factory.set_appearance_rectangle(width=(40, 10, 30), proportion=0.5,
                                         attributes=distr.Levels(["blue", "green"],
                                                                 exact_weighting=True))
        self.rect_stim = factory.create_random_array(n_objects=N)

        try:
            os.mkdir(TEMP_FLD)
        except FileExistsError:
            pass

    def make_path(self, name):
        rtn = os.path.join(TEMP_FLD, name)
        return rtn


class SaveLoad(TestArrays):

    def test_save_load_json(self):
        flname = self.make_path("dots.json")
        self.dot_stim.save(flname)
        new_dots = pynsn.DotArray.load(flname)

        flname = self.make_path("rects.json")
        self.rect_stim.save(flname)
        new_rects = pynsn.RectangleArray.load(flname)

        for flag in [fit.flags.TOTAL_PERIMETER,
                     fit.flags.AV_SURFACE_AREA,
                     fit.flags.SPARSITY,
                     fit.flags.FIELD_AREA,
                     fit.flags.FIELD_AREA_POSITIONS,
                     fit.flags.NUMEROSITY]:

            self.assertAlmostEqual(new_dots.properties.get(flag),
                                   self.dot_stim.properties.get(flag))

            self.assertAlmostEqual(new_rects.properties.get(flag),
                                   self.rect_stim.properties.get(flag))


