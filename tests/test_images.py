from test_save_load import TestArrays
# TODO pygame test missing


class Images(TestArrays):

    def test_pil_images(self):
        from PIL.Image import Image as PILImage
        from pynsn.image import pil_image

        img = pil_image.create(self.dot_stim, self.my_colours)
        self.assertTrue(isinstance(img, PILImage))
        img.save(self.make_path("dots.png"))

        img = pil_image.create(self.rect_stim, self.my_colours)
        self.assertTrue(isinstance(img, PILImage))
        img.save(self.make_path("rects.png"))

    def test_matlab_images(self):
        from matplotlib.pylab import Figure
        from pynsn.image import mpl_figure

        img = mpl_figure.create(self.dot_stim, self.my_colours)
        self.assertTrue(isinstance(img, Figure))
        img.savefig(self.make_path("dots_matplot.png"))

        img = mpl_figure.create(self.rect_stim, self.my_colours)
        self.assertTrue(isinstance(img, Figure))
        img.savefig(self.make_path("rects_matplot.png"))

    def test_SVG_images(self):
        from svgwrite.drawing import Drawing
        from pynsn.image import svg_file

        img = svg_file.create(self.dot_stim, self.my_colours,
                              filename=self.make_path("dots.svg"))
        self.assertTrue(isinstance(img, Drawing))
        img.save()
        img = svg_file.create(self.rect_stim, self.my_colours,
                              filename=self.make_path("rects.svg"))
        self.assertTrue(isinstance(img, Drawing))
        img.save()

