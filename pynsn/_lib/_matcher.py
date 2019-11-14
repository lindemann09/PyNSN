from copy import copy
import random
from hashlib import md5
import json

import numpy as np
from scipy import spatial

from . import _misc, _geometry
from ._convex_hull import EfficientConvexHullDots
from ._item_attributes import ItemAttributes
from ._dot_array_features import DotArrayFeatures
from ._shape import Dot
from . import features
from ._dot_array import DotArray

class ArrayFeatureMatcher(object):

    def __init__(self, dot_array):
        self.da = dot_array


    def match(self, match_feature, match_dot_array=None):
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

        assert isinstance(match_feature, features.ALL_VISUAL_FEATURES)

        features.check_feature_list([match_feature],
                                    check_set_value=match_dot_array
                                                             is None)

        # copy and change values to match this stimulus
        match_feat = copy(match_feature)
        if match_dot_array is not None:
            match_feat.adapt_value(match_dot_array)

        # Adapt
        if isinstance(match_feat, features.ItemDiameter):
            self._match_item_diameter(mean_item_diameter=match_feat.value)

        elif isinstance(match_feat, features.ItemPerimeter):
            self._match_item_diameter(mean_item_diameter=match_feat.value/np.pi)

        elif isinstance(match_feat, features.TotalPerimeter):
            mean_dot_diameter = match_feat.value / (self.da.feature.numerosity * np.pi)
            self._match_item_diameter(mean_dot_diameter)

        elif isinstance(match_feat, features.ItemSurfaceArea):
            ta = self.da.feature.numerosity * match_feat.value
            self._match_total_surface_area(surface_area=ta)

        elif isinstance(match_feat, features.TotalSurfaceArea):
            self._match_total_surface_area(surface_area=match_feat.value)

        elif isinstance(match_feat, features.LogSize):
            logtsa = 0.5 * match_feat.value + 0.5 * _misc.log2(self.da.feature.numerosity)
            self._match_total_surface_area(2 ** logtsa)

        elif isinstance(match_feat, features.LogSpacing):
            logfa = 0.5 * match_feat.value + 0.5 * _misc.log2(self.da.feature.numerosity)
            self._match_field_area(field_area=2 ** logfa,
                                   precision=match_feat.spacing_precision)

        elif isinstance(match_feat, features.Sparsity):
            fa = match_feat.value * self.da.feature.numerosity
            self._match_field_area(field_area=fa,
                                   precision=match_feat.spacing_precision)

        elif isinstance(match_feat, features.FieldArea):
            self._match_field_area(field_area=match_feat.value,
                                   precision=match_feat.spacing_precision)

        elif isinstance(match_feat, features.Coverage):
            self._match_coverage(coverage=match_feat.value,
                                 precision=match_feat.spacing_precision,
                                 match_FA2TA_ratio=match_feat.match_ratio_fieldarea2totalarea) #FIXME experimemtal



    def _match_total_surface_area(self, surface_area):
        # changes diameter
        a_scale = (surface_area / self.da.feature.total_surface_area)
        self.da._diameters = np.sqrt(self.da.surface_areas * a_scale) * 2 / np.sqrt(
            np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
        self.da.set_array_modified()

    def _match_item_diameter(self, mean_item_diameter):
        # changes diameter

        scale = mean_item_diameter / self.da.feature.mean_item_diameter
        self.da._diameters = self.da._diameters * scale
        self.da.set_array_modified()


    # some parameter for matching field arrea
    _ITERATIVE_CONVEX_HULL_MODIFICATION = False  # matching convexhull
    _TAKE_RANDOM_DOT_FROM_CONVEXHULL = False  # todo needs testing
    def _match_field_area(self, field_area,
                          precision=features._DEFAULT_SPACING_PRECISION,
                          use_scaling_only=False):
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

        if self.da.feature.field_area is None:
            return  # not defined
        elif field_area > self.da.feature.field_area or use_scaling_only:
            # field area is currently too small or scaling is enforced
            return self.__scale_field_area(field_area=field_area,
                                           precision=precision)
        elif field_area < self.da.feature.field_area:
            # field area is too large
            self.__decrease_field_area_by_replacement(
                    max_field_area=field_area,
                    iterative_convex_hull_modification=
                    ArrayFeatureMatcher._ITERATIVE_CONVEX_HULL_MODIFICATION)
            # ..and rescaling to avoid to compensate for possible too
            # strong decrease
            return self.__scale_field_area(field_area=field_area,
                                           precision=precision)
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
            while self.da.feature.field_area > max_field_area:
                # remove one random outer dot and remember it
                indices = self.da.convex_hull.indices
                if not ArrayFeatureMatcher._TAKE_RANDOM_DOT_FROM_CONVEXHULL:
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

    def __scale_field_area(self, field_area,
                           precision=features._DEFAULT_SPACING_PRECISION):
        """change the convex hull area to a desired size by scale the polar
        positions  with certain precision

        iterative method can takes some time.

        Note: see doc string `_match_field_area`
        """

        current = self.da.feature.field_area

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
            current = self.da.feature.field_area

            if (current < field_area and step < 0) or \
                    (current > field_area and step > 0):
                step *= -0.2  # change direction and finer grain

        self.da._xy = self.da._xy + old_center
        self.da.set_array_modified()

    def _match_coverage(self, coverage, precision=features._DEFAULT_SPACING_PRECISION,
                        match_FA2TA_ratio=0.5):

        # FIXME check drifting outwards if extra space is small and match_FA2TA_ratio=1
        # FIXME when to realign, realignment changes field_area!
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """

        print("WARNING: _match_coverage is a experimental ")
        # dens = convex_hull_area / total_surface_area
        if match_FA2TA_ratio < 0 or match_FA2TA_ratio > 1:
            match_FA2TA_ratio = 0.5

        total_area_change100 = (coverage * self.da.feature.field_area) - self.da.feature.total_surface_area
        d_change_total_area = total_area_change100 * (1 - match_FA2TA_ratio)
        if abs(d_change_total_area) > 0:
            self._match_total_surface_area(surface_area=self.da.feature.total_surface_area + d_change_total_area)

        self._match_field_area(field_area=self.da.feature.total_surface_area / coverage,
                               precision=precision)

