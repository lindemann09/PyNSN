import random
import numpy as np
from ..lib import misc, geometry
from .visual_features import VisualFeatures

class FeatureMatcher(object):
    ITERATIVE_CONVEX_HULL_MODIFICATION = False
    TAKE_RANDOM_DOT_FROM_CONVEXHULL = True  # TODO maybe method for modification, set via GUI
    DEFAULT_SPACING_PRECISION = 0.0001
    DEFAULT_MATCH_FA2TA_RATIO = 0.5

    def __init__(self, dot_array):
        self.da = dot_array

    def set(self, iterative_convex_hull_modification=None,
            take_random_dot_from_convexhull=None,
            default_spacing_precision=None,
            default_match_fa2ta_ratio=None):
        """Changing class settings of feature matcher.

        This changes the settings of all feature matcher.

        :param default_match_fa2ta_ratio:
        :param default_spacing_precision:
        :param iterative_convex_hull_modification:
                matching convex hull TODO docu
        :param take_random_dot_from_convexhull:
                matching log spacing TODO docu
        """
        if isinstance(iterative_convex_hull_modification, bool):
            FeatureMatcher.ITERATIVE_CONVEX_HULL_MODIFICATION = iterative_convex_hull_modification
        if isinstance(take_random_dot_from_convexhull, bool):
            FeatureMatcher.TAKE_RANDOM_DOT_FROM_CONVEXHULL = take_random_dot_from_convexhull
        if isinstance(default_spacing_precision, float):
            FeatureMatcher.DEFAULT_SPACING_PRECISION = default_spacing_precision
        if isinstance(default_match_fa2ta_ratio, float):
            FeatureMatcher.DEFAULT_MATCH_FA2TA_RATIO = default_match_fa2ta_ratio


    def item_diameter(self, value):
        # changes diameter

        scale = value / self.da._features.mean_item_diameter
        self.da._diameters = self.da._diameters * scale
        self.da._features.reset()

    def total_surface_area(self, value):
        # changes diameter
        a_scale = (value / self.da._features.total_surface_area)
        self.da._diameters = np.sqrt(
            self.da.surface_areas * a_scale) * 2 / np.sqrt(
            np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
        self.da._features.reset()

    def field_area(self, value, precision=None, use_scaling_only=False):
        """changes the convex hull area to a desired size with certain precision

        uses scaling radial positions if field area has to be increased
        uses replacement of outer points (and later re-scaling)

        iterative method can takes some time.
        """

        # PROCEDURE
        #
        # increasing field area:
        #   * merely scales polar coordinates of Dots
        #   * uses __scale_field_area
        #
        # decreasing field area:
        #   a) iterative convex hull modification
        #      1. iteratively replacing outer dots to the side (random  pos.)
        #         (resulting FA is likely to small)
        #      2. increase FA by scaling to match precisely
        #         inside the field area
        #      - this methods results in very angular dot arrays, because it
        #           prefers a solution with a small number of convex hull
        #           dots
        #   b) eccentricity criterion
        #      1. determining circle with the required field area
        #      2. replacing all dots outside this circle to the inside
        #         (random pos.) (resulting FA is likely to small)
        #      3. increase FA by scaling to match precisely
        #      - this method will result is rather circular areas

        if precision is None:
            precision = FeatureMatcher.DEFAULT_SPACING_PRECISION

        if self.da._features.field_area is None:
            return  # not defined
        elif value > self.da._features.field_area or use_scaling_only:
            # field area is currently too small or scaling is enforced
            return self.__scale_field_area(field_area=value,
                                           precision=precision)
        elif value < self.da._features.field_area:
            # field area is too large
            self.__decrease_field_area_by_replacement(
                max_field_area=value,
                iterative_convex_hull_modification=
                FeatureMatcher.ITERATIVE_CONVEX_HULL_MODIFICATION)
            # ..and rescaling to avoid to compensate for possible too
            # strong decrease
            return self.__scale_field_area(field_area=value, precision=precision)
        else:
            return

    def __decrease_field_area_by_replacement(self, max_field_area,
                                             iterative_convex_hull_modification):
        """decreases filed area by recursively moving the most outer point
        to some more central free position (avoids overlapping)

        return False if not possible else True

        Note: see doc string `_match_field_area`

        """

        # centered points
        old_center = self.da.center_of_outer_positions
        self.da._xy = self.da._xy - old_center

        removed_dots = []

        if iterative_convex_hull_modification:
            while self.da._features.field_area > max_field_area:
                # remove one random outer dot and remember it
                indices = self.da._features.convex_hull.indices
                if not FeatureMatcher.TAKE_RANDOM_DOT_FROM_CONVEXHULL:
                    # most outer dot from convex hull
                    radii_outer_dots = geometry.cartesian2polar(self.da._xy[indices],
                                                                radii_only=True)
                    i = np.where(radii_outer_dots == max(radii_outer_dots))[0]
                    idx = indices[i][0]
                else:
                    # remove random
                    idx = indices[random.randint(0, len(indices) - 1)]

                removed_dots.extend(self.da.get(indices=[idx]))
                self.da.delete(idx)

            # add dots to free pos inside the convex hall
            for d in removed_dots:
                try:
                    d._xy = self.da.random_free_position(d.diameter,
                                                         allow_overlapping=False,
                                                         prefer_inside_field_area=True)
                except:
                    raise RuntimeError("Can't find a free position while decreasing field area.\n" + \
                                       "n={}; current FA={}, max_FA={}".format(
                                           self.da._features.numerosity + 1,
                                           self.da._features.field_area,
                                           max_field_area))

                self.da.append_dot(d)

        else:
            # eccentricity criterion
            max_radius = np.sqrt(max_field_area / np.pi)  # for circle with
            # required FA
            idx = np.where(geometry.cartesian2polar(self.da._xy, radii_only=True) > max_radius)[0]
            removed_dots.extend(self.da.get(indices=idx))
            self.da.delete(idx)

            # add inside the circle
            min_dist = self.da.target_array_radius - max_radius + 1
            for d in removed_dots:
                try:
                    d._xy = self.da.random_free_position(d.diameter,
                                                         prefer_inside_field_area=False,
                                                         allow_overlapping=False,
                                                         min_distance_area_boarder=min_dist)
                except:
                    raise RuntimeError(
                        "Can't find a free position while decreasing field area.\n" + \
                        "n={}; current FA={}, max_FA={}".format(
                            self.da._features.numerosity + 1,
                            self.da._features.field_area,
                            max_field_area))
                self.da.append_dot(d)

        self.da._xy = self.da._xy + old_center
        self.da._features.reset()

    def __scale_field_area(self, field_area, precision):
        """change the convex hull area to a desired size by scale the polar
        positions  with certain precision

        iterative method can takes some time.

        Note: see doc string `_match_field_area`
        """

        current = self.da._features.field_area

        if current is None:
            return  # not defined

        scale = 1  # find good scale
        step = 0.1
        if field_area < current:  # current too larger
            step *= -1

        # centered points
        old_center = self.da.center_of_outer_positions
        self.da._xy = self.da._xy - old_center
        centered_polar = geometry.cartesian2polar(self.da._xy)

        # iteratively determine scale
        while abs(current - field_area) > precision:

            scale += step

            self.da._xy = geometry.polar2cartesian(centered_polar * [scale, 1])
            self.da._features.reset()  # required to recalc convex hull
            current = self.da._features.field_area

            if (current < field_area and step < 0) or \
                    (current > field_area and step > 0):
                step *= -0.2  # change direction and finer grain

        self.da._xy = self.da._xy + old_center
        self.da._features.reset()

    def coverage(self, value,
                 precision=None,
                 match_FA2TA_ratio=None):

        # FIXME check drifting outwards if extra space is small and match_FA2TA_ratio=1
        # FIXME when to realign, realignment changes field_area!
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """

        print("WARNING: _match_coverage is a experimental ")
        # dens = convex_hull_area / total_surface_area
        if match_FA2TA_ratio is None:
            match_FA2TA_ratio = FeatureMatcher.DEFAULT_MATCH_FA2TA_RATIO
        elif match_FA2TA_ratio < 0 or match_FA2TA_ratio > 1:
            match_FA2TA_ratio = 0.5
        if precision is None:
            precision = FeatureMatcher.DEFAULT_SPACING_PRECISION

        total_area_change100 = (value * self.da._features.field_area) - self.da._features.total_surface_area
        d_change_total_area = total_area_change100 * (1 - match_FA2TA_ratio)
        if abs(d_change_total_area) > 0:
            self.total_surface_area(self.da._features.total_surface_area + d_change_total_area)

        self.field_area(self.da._features.total_surface_area / value,
                        precision=precision)

    def item_perimeter(self, value):

        self.item_diameter(value / np.pi)

    def total_perimeter(self, value):
        tmp = value / (self.da._features.numerosity * np.pi)
        self.item_diameter(tmp)

    def item_surface_area(self, value):
        ta = self.da._features.numerosity * value
        self.total_surface_area(ta)

    def log_spacing(self, value, precision=None):

        logfa = 0.5 * value + 0.5 * misc.log2(
            self.da._features.numerosity)
        self.field_area(value=2 ** logfa, precision=precision)

    def log_size(self, value):
        logtsa = 0.5 * value + 0.5 * misc.log2(self.da._features.numerosity)
        self.total_surface_area(2 ** logtsa)

    def sparcity(self, value, precision=None):
        self.field_area(value=value * self.da._features.numerosity,
                        precision=precision)

    def match_feature(self, feature, value):
        """
        match_properties: continuous property or list of continuous properties
        several properties to be matched
        if match dot array is specified, array will be match to match_dot_array, otherwise
        the values defined in match_features is used.
        some matching requires realignement to avoid overlaps. However,
        realigment might result in a different field area. Thus, realign after
        matching for  Size parameter and realign before matching space
        parameter.

        """

        # type check
        if feature not in VisualFeatures.ALL_FEATURES:
            raise ValueError("{} is not a visual feature.".format(feature))

        # Adapt
        if feature == VisualFeatures.ITEM_DIAMETER:
            self.item_diameter(value=value)

        elif feature == VisualFeatures.ITEM_PERIMETER:
            self.item_perimeter(value=value)

        elif feature == VisualFeatures.TOTAL_PERIMETER:
            self.total_perimeter(value=value)

        elif feature == VisualFeatures.ITEM_SURFACE_AREA:
            self.item_surface_area(value=value)

        elif feature == VisualFeatures.TOTAL_SURFACE_AREA:
            self.total_surface_area(value=value)

        elif feature == VisualFeatures.LOG_SIZE:
            self.log_size(value=value)

        elif feature == VisualFeatures.LOG_SPACING:
            self.log_spacing(value=value)

        elif feature == VisualFeatures.SPARSITY:
            self.sparcity(value=value)

        elif feature == VisualFeatures.FIELD_AREA:
            self.field_area(value=value)

        elif feature == VisualFeatures.COVERAGE:
            self.coverage(value=value)
