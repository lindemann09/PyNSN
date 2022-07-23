from os import path
import unittest
import tempfile

from pynsn import NSNFactory, ImageColours
from pynsn import distributions as distr

TEMP_FLD = tempfile.mkdtemp(prefix="pynsn_images_")
print("Image folder {}".format(TEMP_FLD))


class Images(unittest.TestCase):

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
                              ) # FIXME test opacity

        factory = NSNFactory(target_area_radius=200)
        factory.set_appearance_dot(diameter=(40, 10, 30),
                                        attributes=distr.Levels(["blue", "green"],
                                        exact_weighting=True))
        self.dot_stim = factory.create_random_array(n_objects=N)

        factory.set_appearance_rectangle(width=(40, 10, 30), proportion=0.5,
                                                attributes=distr.Levels(["blue", "green"],
                                                exact_weighting=True))
        self.rect_stim = factory.create_random_array(n_objects=N)

    def make_path(self, name):
        rtn = path.join(TEMP_FLD, name)
        return rtn

    def test_pil_images(self):
        from PIL.Image import Image as PILImage
        from pynsn.image import pil_image

        img = pil_image.create(self.dot_stim, self.my_colours)
        assert isinstance(img, PILImage)
        img.save(self.make_path("dots.png"))

        img = pil_image.create(self.rect_stim, self.my_colours)
        assert isinstance(img, PILImage)
        img.save(self.make_path("rects.png"))

    def test_matlab_images(self):
        from matplotlib.pylab import Figure
        from pynsn.image import mpl_figure

        img = mpl_figure.create(self.dot_stim, self.my_colours)
        assert isinstance(img, Figure)
        img.savefig(self.make_path("dots_matplot.png"))

        img = mpl_figure.create(self.rect_stim, self.my_colours)
        assert isinstance(img, Figure)
        img.savefig(self.make_path("rects_matplot.png"))

    def test_SVG_images(self):
        from svgwrite.drawing import Drawing
        from pynsn.image import svg_file

        img = svg_file.create(self.dot_stim, self.my_colours,
                              filename=self.make_path("dots.svg"))
        assert isinstance(img, Drawing)
        img.save()

        img = svg_file.create(self.rect_stim, self.my_colours,
                              filename=self.make_path("rects.svg"))
        assert isinstance(img, Drawing)
        img.save()
