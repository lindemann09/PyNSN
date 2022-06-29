import unittest
from pynsn import arrays, distr, random_array, match, VisualFeature
from test_dot_array import DotsSmall


class RectanglesSmall(DotsSmall):
    def settings(self):
        super().settings()
        self.size_dist = random_array.SizeDistribution(
            width=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
            height=distr.Normal(min_max=(10, 40), mu=20, sigma=10))

# class RectanglesMedium(DotsMedium):
#     def settings(self):
#         super().settings()
#         self.size_dist = random_array.SizeDistribution(
#             width=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
#             height=distr.Normal(min_max=(10, 40), mu=20, sigma=10))
#
# class RectanglesLarge(DotsLarge):
#     def settings(self):
#         super().settings()
#         self.size_dist = random_array.SizeDistribution(
#             width=distr.Normal(min_max=(5, 30), mu=10, sigma=5),
#             height=distr.Normal(min_max=(5, 30), mu=10, sigma=5))


if __name__ == "__main__":
    unittest.main()
