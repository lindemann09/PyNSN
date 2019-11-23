from copy import copy
import random
import numpy as np
from . import _misc, _geometry
from ._visual_features import Features

DEFAULT_SPACING_PRECISION = 0.0001
DEFAULT_MATCH_FA2TA_RATIO = 0.5

# some parameter for matching field arrea
ITERATIVE_CONVEX_HULL_MODIFICATION = False  # matching convexhull TODO DOCU
# and matching log spacing
TAKE_RANDOM_DOT_FROM_CONVEXHULL = False  # todo needs testing


class FeatureMatcher(object):

    def __init__(self, dot_array):
        self.da = dot_array

    def item_diameter(self, value):
        # changes diameter

        scale = value / self.da.features.mean_item_diameter
        self.da._diameters = self.da._diameters * scale
        self.da.set_array_modified()

    def total_surface_area(self, value):
        # changes diameter
        a_scale = (value / self.da.features.total_surface_area)
        self.da._diameters = np.sqrt(
            self.da.surface_areas * a_scale) * 2 / np.sqrt(
            np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
        self.da.set_array_modified()


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
        #      - this method will result is rather circulr areas


        if precision is None:
            precision = DEFAULT_SPACING_PRECISION

        if self.da.features.field_area is None:
            return  # not defined
        elif value > self.da.features.field_area or use_scaling_only:
            # field area is currently too small or scaling is enforced
            return self.__scale_field_area(field_area=value,
                                           precision=precision)
        elif value < self.da.features.field_area:
            # field area is too large
            self.__decrease_field_area_by_replacement(
                    max_field_area=value,
                    iterative_convex_hull_modification=
                    ITERATIVE_CONVEX_HULL_MODIFICATION)
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
            while self.da.features.field_area > max_field_area:
                # remove one random outer dot and remember it
                indices = self.da.convex_hull.indices
                if not TAKE_RANDOM_DOT_FROM_CONVEXHULL:
                    # most outer dot from convex hull
                    radii_outer_dots = _geometry.cartesian2polar(self.da.xy[indices],
                                                             radii_only=True)
                    i = np.where(radii_outer_dots==max(radii_outer_dots))[0]
                    idx = indices[i][0]
                else:
                    # remove random
                    idx = indices[random.randint(0, len(indices)-1)]

                removed_dots.extend(self.da.get_dots(indices=[idx]))
                self.da.delete(idx)

            # add dots to free pos inside the convex hall
            for d in removed_dots:
                d.xy = self.da.random_free_dot_position(d.diameter,
                                                   allow_overlapping=False,
                                                   prefer_inside_field_area=True)
                self.da.append_dot(d)

        else:
            # eccentricity criterion
            max_radius =  np.sqrt(max_field_area/np.pi) # for circle with
                                                        # required FA
            idx = np.where(_geometry.cartesian2polar(self.da.xy, radii_only=True) > max_radius)[0]
            removed_dots.extend(self.da.get_dots(indices=idx))
            self.da.delete(idx)

            # add inside the circle
            min_dist = self.da.target_array_radius - max_radius + 1
            for d in removed_dots:
                d.xy = self.da.random_free_dot_position(d.diameter,
                                            allow_overlapping=False,
                                            min_distance_area_boarder=min_dist)
                self.da.append_dot(d)

        self.da._xy = self.da._xy + old_center
        self.da.set_array_modified()

    def __scale_field_area(self, field_area, precision):
        """change the convex hull area to a desired size by scale the polar
        positions  with certain precision

        iterative method can takes some time.

        Note: see doc string `_match_field_area`
        """

        current = self.da.features.field_area

        if current is None:
            return  # not defined

        scale = 1  # find good scale
        step = 0.1
        if field_area < current:  # current too larger
            step *= -1

        # centered points
        old_center = self.da.center_of_outer_positions
        self.da._xy = self.da._xy - old_center
        centered_polar = _geometry.cartesian2polar(self.da._xy)

        # iteratively determine scale
        while abs(current - field_area) > precision:

            scale += step

            self.da._xy = _geometry.polar2cartesian(centered_polar * [scale, 1])
            self.da.set_array_modified()  # required to recalc convex hull
            current = self.da.features.field_area

            if (current < field_area and step < 0) or \
                    (current > field_area and step > 0):
                step *= -0.2  # change direction and finer grain

        self.da._xy = self.da._xy + old_center
        self.da.set_array_modified()

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
            match_FA2TA_ratio = DEFAULT_MATCH_FA2TA_RATIO
        elif match_FA2TA_ratio < 0 or match_FA2TA_ratio > 1:
            match_FA2TA_ratio = 0.5

        total_area_change100 = (value * self.da.features.field_area) - self.da.features.total_surface_area
        d_change_total_area = total_area_change100 * (1 - match_FA2TA_ratio)
        if abs(d_change_total_area) > 0:
            self.total_surface_area(self.da.features.total_surface_area + d_change_total_area)

        self.field_area(self.da.features.total_surface_area / value,
                        precision=precision)


    def item_perimeter(self, value):

        self.item_diameter(value / np.pi)

    def total_perimeter(self, value):
        tmp = value / (self.da.features.numerosity * np.pi)
        self.item_diameter(tmp)

    def item_surface_area(self, value):
        ta = self.da.features.numerosity * value
        self.total_surface_area(ta)

    def log_spacing(self, value, precision=None):

        logfa = 0.5 * value + 0.5 * _misc.log2(
            self.da.features.numerosity)
        self.field_area(value=2 ** logfa, precision=precision)

    def log_size(self, value):
        logtsa = 0.5 * value + 0.5 * _misc.log2(self.da.features.numerosity)
        self.total_surface_area(2 ** logtsa)

    def sparcity(self, value, precision = None):
        self.field_area(value= value * self.da.features.numerosity,
                        precision=precision)


    def match_feature(self, feature,
                      value = None,
                      reference_dot_array=None):
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

        if value is None and reference_dot_array is None:
            raise ValueError("Please specify a value or a "
                             "reference_dot_array.")

        if value is not None and reference_dot_array is not None:
            raise ValueError("Please specify either a value or "
                             "reference_dot_array, not both.")

        # type check
        if feature not in Features.ALL_FEATURES:
            raise ValueError("{} is not a visual feature.".format(feature))

        if value is None:
            value = reference_dot_array.features.get(feature)

        # Adapt
        if feature == Features.ITEM_DIAMETER:
            self.item_diameter(value=value)

        elif feature == Features.ITEM_PERIMETER:
            self.item_perimeter(value=value)

        elif feature == Features.TOTAL_PERIMETER:
            self.total_perimeter(value=value)

        elif feature == Features.ITEM_SURFACE_AREA:
            self.item_surface_area(value=value)

        elif feature == Features.TOTAL_SURFACE_AREA:
            self.total_surface_area(value=value)

        elif feature == Features.LOG_SIZE:
            self.log_size(value=value)

        elif feature == Features.LOG_SPACING:
            self.log_spacing(value=value,
                             precision=DEFAULT_SPACING_PRECISION)

        elif feature == Features.SPARSITY:
            self.sparcity(value=value, precision=DEFAULT_SPACING_PRECISION)

        elif feature == Features.FIELD_AREA:
            self.field_area(value=value, precision=DEFAULT_SPACING_PRECISION)

        elif feature == Features.COVERAGE:
            self.coverage(value=value, precision=DEFAULT_SPACING_PRECISION)
