import unittest
from pynsn import arrays, distr, random_array, adapt, VisualFeature, scale


class DotsSmall(unittest.TestCase):

    def settings(self):
        self.ref = arrays.GenericObjectArray(target_area_radius=200)
        self.size_dist = random_array.SizeDistribution(
            diameter=distr.Beta(min_max=(10, 30), mu=15, sigma=2))
        self.n_dots = 5

    def setUp(self):
        self.settings()
        self.stimulus = random_array.create(reference_array=self.ref,
                                       size_distribution=self.size_dist,
                                       n_objects=self.n_dots)

    def test_numerosity(self):
        self.assertEqual(self.stimulus.features.numerosity, self.n_dots)

    def change_feature(self, feature, first=0.8, second=1.15,
                       scale_factor = 1.1, places=7):

        # check scale
        stim = self.stimulus.copy()
        new_value = stim.features.get(feature) * scale_factor
        scale.visual_feature(stim, feature, scale_factor)
        self.assertAlmostEqual(stim.features.get(feature),
                               new_value, places=places)

        #  changes two time feature (e.g. first decrease and then increase)
        # first
        new_value = stim.features.get(feature) * first
        adapt.visual_feature(stim, feature, new_value)
        self.assertAlmostEqual(stim.features.get(feature),
                               new_value, places=places)
        #second
        new_value = stim.features.get(feature) * second
        adapt.visual_feature(stim, feature, new_value)
        self.assertAlmostEqual(stim.features.get(feature), new_value,
                               places=places)


    def test_match_av_surface_area(self):
        # decrease
        self.change_feature(feature=VisualFeature.AV_SURFACE_AREA)

    def test_match_av_size(self):
        # decrease
        self.change_feature(feature=VisualFeature.AV_DOT_DIAMETER)

    def test_match_av_perimeter(self):
        # decrease
        self.change_feature(feature=VisualFeature.AV_PERIMETER)

    def test_match_total_surface_area(self):
        # decrease
        self.change_feature(feature=VisualFeature.TOTAL_SURFACE_AREA)

    def test_match_total_perimeter(self):
        # decrease
        self.change_feature(feature=VisualFeature.TOTAL_PERIMETER)

    def test_match_field_area(self):
        # decrease
        self.change_feature(feature=VisualFeature.SPARSITY, first=1.2,
                            second=0.85, places=4)

    def test_match_sparcity(self):
        # decrease
        self.change_feature(feature=VisualFeature.SPARSITY, first=1.2,
                            second=0.85, places=4)

    def test_match_log_size(self):
        # decrease
        self.change_feature(feature=VisualFeature.LOG_SIZE)

    def test_match_log_spacing(self):
        # decrease
        self.change_feature(feature=VisualFeature.LOG_SPACING, first=1.2,
                            second=0.85)


class DotsMedium(DotsSmall):
    def settings(self):
        super().settings()
        self.n_dots = 25

class DotsLarge(DotsSmall):
    def settings(self):
        super().settings()
        self.n_dots = 75
