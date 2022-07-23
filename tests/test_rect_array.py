import unittest
from pynsn import distributions as distr
from pynsn._visual_properties import fit, scale
from test_dot_array import DotsSmall


class RectanglesSmall(DotsSmall):
    def settings(self):
        super().settings()
        self.factory.set_appearance_rectangle(
            width=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
            height=distr.Normal(min_max=(10, 40), mu=20, sigma=10))

    def test_fit_av_size(self):
        first = 0.8
        second = 1.15
        scale_factor = 1.1
        places = 7
        # check scale
        stim = self.stimulus.copy()
        new_value = stim.properties.average_rectangle_size * scale_factor
        scale.average_rectangle_size(stim, scale_factor)
        self.assertAlmostEqual(stim.properties.average_rectangle_size[0],
                               new_value[0], places=places)
        self.assertAlmostEqual(stim.properties.average_rectangle_size[1],
                               new_value[1], places=places)

        #  changes two time feature (e.g. first decrease and then increase)
        # first
        new_value = stim.properties.average_rectangle_size * first
        fit.average_rectangle_size(stim, new_value)
        self.assertAlmostEqual(stim.properties.average_rectangle_size[0],
                               new_value[0], places=places)
        self.assertAlmostEqual(stim.properties.average_rectangle_size[1],
                               new_value[1], places=places)
        #second
        new_value = stim.properties.average_rectangle_size * second
        fit.average_rectangle_size(stim, new_value)
        self.assertAlmostEqual(stim.properties.average_rectangle_size[0],
                               new_value[0], places=places)
        self.assertAlmostEqual(stim.properties.average_rectangle_size[1],
                               new_value[1], places=places)


class RectanglesMedium(RectanglesSmall):
    def settings(self):
        super().settings()
        self.factory.set_appearance_rectangle(
            width=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
            height=distr.Normal(min_max=(10, 40), mu=20, sigma=10))
        self.n_dots = 25


class RectanglesLarge(RectanglesSmall):
    def settings(self):
        super().settings()
        self.factory.set_appearance_rectangle(
            width=distr.Normal(min_max=(5, 30), mu=10, sigma=5),
            height=distr.Normal(min_max=(5, 30), mu=10, sigma=5))
        self.n_dots = 75


if __name__ == "__main__":
    unittest.main()
