"""
Dot Array
"""
from __future__ import print_function, division, unicode_literals

from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import random
from copy import deepcopy, copy
import numpy as np
from .dot import Dot
from .dot_list import DotList, DotListProperties, numpy_vector
from .colours import convert_colour

TWO_PI = 2 * np.pi


class DotArray(DotList):

    def __init__(self, max_array_radius, minimum_gap=1, xy=None,
                 diameters=None, colours=None, pictures=None):
        """Dot array is restricted to a certain area and can generate random dots and
            be realigned """

        DotList.__init__(self)
        self.max_array_radius = max_array_radius
        self.minimum_gap = minimum_gap
        self._colours = np.array([])
        self._pictures = np.array([])
        if (xy, diameters, colours, pictures) != (None, None, None, None):
            self.append(xy, diameters, colours=colours, pictures=pictures)

    @property
    def colours(self):
        return self._colours

    @property
    def pictures(self):
        return self._pictures

    def clear(self):
        DotList.clear(self)
        self._colours = np.array([])
        self._pictures = np.array([])

    def append(self, xy, diameters, colours=None, pictures=None):
        """append dots using numpy array"""

        # ensure good numpy array or None
        diameters = numpy_vector(diameters)
        pictures = numpy_vector(pictures)
        colours = numpy_vector(colours)
        if (colours is not None and len(colours) != len(diameters)) or \
                (pictures is not None and len(pictures) != len(diameters)):
            raise RuntimeError(u"Bad shaped data: " + u"colour and/or picture has not the same length as diameter")

        DotList.append(self, xy, diameters)

        if pictures is None:
            pictures = [None]*len(diameters)
        self._pictures = np.append(self._pictures, pictures)

        if colours is None:
            colours = [None]*len(diameters)
        for c in colours:
            self._colours = np.append(self._colours, convert_colour(c))

    def delete(self, index):
        DotList.delete(self, index)
        self._colours = np.delete(self._colours, index)
        self._pictures = np.delete(self._pictures, index)

    def copy(self, indices=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """
        if indices is None:
            indices = list(range(self.prop_numerosity))

        return DotArray(max_array_radius=self.max_array_radius,
                        minimum_gap=self.minimum_gap,
                        xy=self._xy[indices, :].copy(),
                        diameters=self._diameters[indices].copy(),
                        colours=self._colours[indices].copy(),
                        pictures=self._pictures[indices].copy())

    def append_dot(self, dot):
        self.append(xy=[dot.x, dot.y], diameters=dot.diameter,
                    colours=dot.colour, pictures=dot.picture)

    def join(self, dot_array, realign=True):
        """add another dot arrays"""

        self.append(xy=dot_array._xy, diameters=dot_array._diameters,
                    colours=dot_array._colours, pictures=dot_array._pictures)
        if realign:
            self.realign()

    def get_dots(self, indices=None, diameter=None, colour=None, picture=None):
        """returns all dots
         filtering possible:
         if diameter/colour/picture is defined it returns only dots a particular diameter/colour/picture

         """
        rtn = []
        i = -1
        if indices is not None:
            try:
                indices = list(indices)  # check if iterable
            except:
                indices = [indices]

        for xy, dia, col, pict in zip(self._xy, self._diameters,
                                      self._colours, self._pictures):
            i += 1
            if (indices is not None and i not in indices) or \
                    (diameter is not None and dia != diameter) or \
                    (colour is not None and col != colour) or \
                    (picture is not None and pict != picture):
                continue

            rtn.append(Dot(x=xy[0], y=xy[1], diameter=dia,
                           colour=col, picture=pict))
        return rtn

    def change_colour(self, colour, indices=None):
        """allows usung color names and rgb array, since it
        converts colour """

        if isinstance(indices, int):
            indices = [indices]
        elif indices is None:
            indices = range(len(self._colours))

        self._colours[indices] = convert_colour(colour)

    def realign(self):
        """Realigns the dots in order to remove all dots overlaps. If two dots
        overlap, the dots that is further apart from the arry center will be
        moved opposite to the direction of the other dot until there is no
        overlap (note: minimun_gap parameter). If two dots have exactly the same
        position the same position one is minimally shifted in a random direction.

        """

        shift_required = False
        error = False

        # from inner to outer remove overlaps
        for i in np.argsort(DotList._cartesian2polar(self._xy, radii_only=True)):
            if self._remove_overlap_for_dot(dot_id=i, minimum_gap=self.minimum_gap):
                shift_required = True

        # sqeeze in points that pop out of the stimulus area radius
        cnt = 0
        while True:
            radii = DotList._cartesian2polar(self._xy, radii_only=True)
            too_far = np.where((radii + self._diameters // 2) > self.max_array_radius)[0]  # find outlier
            if len(too_far) > 0:

                # squeeze in outlier
                polar = DotList._cartesian2polar([self._xy[too_far[0], :]])[0]
                polar[0] = self.max_array_radius - self._diameters[
                    too_far[0]] // 2 - 0.000000001  # new radius #todo check if 0.00001 required
                new_xy = DotList._polar2cartesian([polar])[0]
                self._xy[too_far[0], :] = new_xy

                # remove overlaps centered around new outlier position
                self._xy -= new_xy
                # remove all overlaps (inner to outer, i.e. starting with outlier)
                for i in np.argsort(DotList._cartesian2polar(self._xy, radii_only=True)):
                    self._remove_overlap_for_dot(dot_id=i, minimum_gap=self.minimum_gap)
                # new pos for outlyer
                self._xy += new_xy  # move back to old position
                shift_required = True
            else:
                break  # end while loop

            cnt += 1
            if cnt > 20:
                error = True
                break

        if error:
            return False, u"Can't find solution when removing outlier (n=" + \
                   str(self.prop_numerosity) + ")"
        if not shift_required:
            return True, ""
        else:
            return self.realign()  # recursion

    @property
    def prop_density_max_field(self):
        """density takes into account the full possible dot area """
        return np.pi * self.max_array_radius ** 2 / self.prop_total_surface_area

    def split_array_by_colour(self):
        """returns a list of arrays
        each array contains all dots of with particular colour"""
        rtn = []
        for c in np.unique(self._colours):
            if c is not None:
                da = DotArray(max_array_radius=self.max_array_radius,
                              minimum_gap=self.minimum_gap)
                for d in self.get_dots(colour=c):
                    da.append_dot(d)
                rtn.append(da)
        return rtn

    def get_properties_split_by_colours(self):
        """returns is unicolor or no color"""
        if len(np.unique(self._colours)) == 1:
            return None

        rtn = DotListProperties._make_arrays()
        for da in self.split_array_by_colour():
            prop = da.get_properties()
            rtn[0].append(self.object_id + str(da._colours[0]))
            for i in range(1, len(rtn)):
                rtn[i].append(prop[i])
        return rtn

    def __str__(self):
        return self.get_csv()

    def get_csv(self, num_format="%7.2f", variable_names=True,
                object_id_column=True, num_idx_column=True,
                colour_column=False, picture_column=False):
        """Return the dot array as csv text

        Parameter
        ---------
        variable_names : bool, optional
            if True variable name will be printed in the first line

        """

        rtn = ""
        if variable_names:
            if object_id_column:
                rtn += u"object_id,"
            if num_idx_column:
                rtn += u"num_id,"
            rtn += u"x,y,diameter"
            if colour_column:
                rtn += u",colour"
            if picture_column:
                rtn += u",file"
            rtn += u"\n"

        if object_id_column:
            obj_id = self.object_id

        for cnt in range(len(self._xy)):
            if object_id_column:
                rtn += "{0}, ".format(obj_id)
            if num_idx_column:
                rtn += "{},".format(self.prop_numerosity)
            rtn += num_format % self._xy[cnt, 0] + "," + num_format % self._xy[cnt, 1] + "," + \
                   num_format % self._diameters[cnt]
            if colour_column:
                rtn += ", {}".format(self._colours[cnt])
            if picture_column:
                rtn += ", {}".format(self._pictures[cnt])
            rtn += "\n"
        return rtn

    def random_free_dot_position(self, dot_diameter,
                                 ignore_overlapping=False):
        """moves a dot to an available random position

        raise exception if not found"""

        cnt = 0
        while True:
            cnt += 1
            proposal_polar = np.random.random(2) * \
                             (self.max_array_radius - (dot_diameter / 2.0), TWO_PI)
            proposal_xy = DotList._polar2cartesian([proposal_polar])[0]

            bad_position = False
            if not ignore_overlapping:
                # find bad_positions
                dist = self.distances(proposal_xy, dot_diameter)
                idx = np.where(dist < self.minimum_gap)[0]  # overlapping dot ids
                bad_position = len(idx) > 0

            if not bad_position:
                return proposal_xy
            elif cnt > 3000:
                raise RuntimeError(u"Can't find a solution")

    def shuffle_all_positions(self, ignore_overlapping=False):
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        diameters = self._diameters
        self._diameters = np.array([])
        self._xy = np.array([])
        for d in diameters:
            xy = self.random_free_dot_position(d, ignore_overlapping=ignore_overlapping)
            self._diameters = np.append(self._diameters, d)
            if len(self._xy) == 0:
                self._xy = np.array([xy])
            else:
                self._xy = np.append(self._xy, [xy], axis=0)

    def number_deviant(self, change_numerosity):
        """number deviant
        """

        # make a copy for the deviant
        deviant = self.copy()
        if self.prop_numerosity + change_numerosity <= 0:
            deviant.clear()
        else:
            # add or remove random dots
            for _ in range(abs(change_numerosity)):
                rnd = random.randint(0,
                                     deviant.prop_numerosity - 1)  # strange: do not use np.rand, because parallel running process produce the same numbers
                if change_numerosity < 0:
                    # remove dots
                    deviant.delete(rnd)
                else:
                    # copy a random dot
                    deviant.append(xy=self.random_free_dot_position(dot_diameter=self._diameters[rnd]),
                                   diameters=self._diameters[rnd],
                                   colours=self._colours[rnd],
                                   pictures=self._pictures[rnd])

        return deviant

    def match_total_surface_area(self, surface_area):
        # changes diameter
        a_scale = (surface_area / self.prop_total_surface_area)
        self._diameters = np.sqrt(self.surface_areas * a_scale) * 2 / np.sqrt(
            np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)

    def match_mean_dot_diameter(self, mean_dot_diameter):
        # changes diameter

        scale = mean_dot_diameter / self.prop_mean_dot_diameter
        self._diameters = self._diameters * scale

    def match_total_circumference(self, total_circumference):
        # linear to fit_mean_dot_diameter, but depends on numerosity
        mean_dot_diameter = total_circumference / (self.prop_numerosity * np.pi)
        self.match_mean_dot_diameter(mean_dot_diameter)

    def match_convex_hull_area(self, convex_hull_area, precision=0.001,
                               center_array=True):
        """changes the convex hull area to a desired size with certain precision

        iterative method can takes some time.
        """

        current = self.prop_convex_hull_area

        if current is None:
            return  # not defined

        # iteratively determine scale
        scale = 1  # find good scale
        step = 0.1
        if convex_hull_area < current:  # current too larger
            step *= -1

        # centered  points
        old_center = self.center_of_outer_positions
        centered_polar = DotList._cartesian2polar(self._xy - old_center)

        while abs(current - convex_hull_area) > precision:
            scale += step

            self._xy = DotList._polar2cartesian(centered_polar * [scale, 1])
            current = self.prop_convex_hull_area

            if (current < convex_hull_area and step < 0) or \
                    (current > convex_hull_area and step > 0):
                step *= -0.2  # change direction and finer grain

        if not center_array:
            self._xy += old_center

    def match_density(self, density, precision=0.001, ratio_convex_hull2area_adaptation=0.5):
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """

        # dens = convex_hull_area / total_surface_area
        if ratio_convex_hull2area_adaptation < 0 or ratio_convex_hull2area_adaptation > 1:
            ratio_convex_hull2area_adaptation = 0.5

        area_change100 = (self.prop_convex_hull_area / density) - self.prop_total_surface_area
        d_change_area = area_change100 * (1 - ratio_convex_hull2area_adaptation)
        if abs(d_change_area) > 0:
            self.match_total_surface_area(surface_area=self.prop_total_surface_area + d_change_area)

        self.match_convex_hull_area(convex_hull_area=self.prop_total_surface_area * density,
                                    precision=precision)
