"""
Dot Array
"""
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from os import listdir
import types
import random
from copy import deepcopy, copy
import numpy as np
from scipy.spatial import ConvexHull
from . import Dot, random_beta, shape_parameter_beta
TWO_PI = 2*np.pi

PROPERTY_FILE_SUFFIX = ".property.csv"

def list_all_saved_incremental_arrays(picture_folder):
    rtn = []
    for x in listdir(picture_folder):
        if x.endswith(PROPERTY_FILE_SUFFIX):
            rtn.append(x[:len(x)-len(PROPERTY_FILE_SUFFIX)])
    return rtn


class NumpyPositionList(object):
    """Numpy Position list

    for optimized calculations
    """

    def __init__(self, list_of_Dots):
        xy = []
        d = []
        for dot in list_of_Dots:
            xy.append(dot.xy)
            d.append(dot.diameter)

        self._xy = np.array(xy)
        self.diameter = np.array(d)
        self._polar = None

    @property
    def xy(self):
        """numpy array of [x, y]"""

        if self._xy is None:
            self.calc_cartesian()
        return self._xy

    @xy.setter
    def xy(self, xy_position_list):
        self._xy = np.array(xy_position_list)
        self._polar = None

    @property
    def polar(self):
        """numpy array of [radius, angle]"""
        if self._polar is None:
            self.calc_polar()
        return self._polar

    @polar.setter
    def polar(self, radius_angle_list):
        self._polar = np.array(radius_angle_list)
        self._xy = None

    def calc_cartesian(self):
        self._xy = np.array([self._polar[:, 0] * np.cos(self._polar[:, 1]),
                             self._polar[:, 0] * np.sin(self._polar[:, 1])]).T

    def calc_polar(self):
        self._polar = np.array([np.hypot(self._xy[:, 0], self._xy[:, 1]),
                                np.arctan2(self._xy[:, 0], self._xy[:, 1])]).T

    def distance(self, dot):
        return np.hypot(self.xy[:,0] - dot.x, self.xy[:, 1] - dot.y) - \
               ((self.diameter + dot.diameter) / 2.0)


    def jitter_identical_positions(self):
        """jitters points with identical position"""
        tmp = Dot(diameter=0)
        for i in range(len(self._xy)):
            ref_dot = Dot(x=self._xy[i, 0], y=self._xy[i, 1], diameter=self.diameter[i])
            identical = np.where(np.all(np.equal(self._xy, ref_dot.xy), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != i:
                        tmp.polar = (0.1, random.random() * TWO_PI)
                        self._xy[x, :] -= tmp.xy
                        self._polar = None

    def remove_overlap_for_dot(self, dotid, min_gap):
        """remove overlap for one point"""

        ref_dot = Dot(x=self._xy[dotid, 0], y=self._xy[dotid, 1], diameter=self.diameter[dotid])
        dist = self.distance(ref_dot)
        tmp = Dot(diameter=0)
        shift_required = False
        for x in np.where(dist < min_gap)[0]: # x=overlapping dot id
            if x != dotid:
                shift_required = True
                # calc vector (tmp) to shift
                tmp.xy = self._xy[x, :] - ref_dot.xy
                tmp.pos_radius = 0.000000001 + min_gap - dist[x]

                self._xy[x, :] = (self._xy[x, 0] + tmp.x, self._xy[x, 1] + tmp.y)
                shift_required = True

        if shift_required:
            self._polar = None

        return shift_required


class DotArray(object):

    def __init__(self, n_dots, stimulus_area_radius,
                 dot_diameter_mean,
                 dot_diameter_range=None,
                 dot_diameter_std=None,
                 dot_picture = None,
                 dot_colour=None,
                 minium_gap=1):

        """Create a Random Dot Kinematogram

        Parameters:
        -----------
        stimulus_area_radius : int
            the radius of the stimulus area
        n_dots : int
            number of moving dots

        FIXME

        Notes:
        ------
        Logging is switch off per default

        """

        if dot_diameter_std <= 0:
            dot_diameter_std = None

        self._min_gap = minium_gap
        self._stimulus_area_radius = stimulus_area_radius
        self._dot_diameter_range = dot_diameter_range
        self._dot_diameter_mean = dot_diameter_mean
        self._dot_diameter_std = dot_diameter_std
        self._dot_colour = dot_colour
        self._dot_picture = dot_picture
        self.dots = []
        self._create_dots(n_dots=n_dots)

    def copy(self, subset_dot_ids=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if subset_dot_ids is None:
            return deepcopy(self)
        else:
            rtn = copy(self)
            rtn.dots = []
            if isinstance(subset_dot_ids, types.IntType):
                subset_dot_ids = [subset_dot_ids]
            for x in subset_dot_ids:
                rtn.dots.append(deepcopy(self.dots[x]))
            return rtn

    def _create_dots(self, n_dots):
        self.dots = []
        for _ in range(n_dots):
            self.dots.append(self.create_random_dot(ignore_overlapping=False))

    def create_random_dot(self, ignore_overlapping=False):
        """returns a random dot

        if no available space can be found it raise an Runtime Error
        """

        if self._dot_diameter_range is None or \
                        self._dot_diameter_std is None:
            # constant mean
            rtn = Dot(diameter=self._dot_diameter_mean,
                              colour=self._dot_colour,
                              picture=self._dot_picture)
        else:
            # draw diameter from beta distribution
            parameter = shape_parameter_beta(self._dot_diameter_range,
                                             self._dot_diameter_mean,
                                             self._dot_diameter_std)
            rtn = Dot(diameter=random_beta(
                self._dot_diameter_range, parameter),
                colour=self._dot_colour)

        self.randomize_dot_position(rtn, ignore_overlapping)
        return rtn

    def randomize_dot_position(self, dot,
                               ignore_overlapping=False,
                               minimum_gap=None):
        """moves a dot to an available random position

        raise exception if not found"""

        if minimum_gap is not None:
            self._min_gap = minimum_gap

        cnt = 0
        while True:
            cnt += 1
            dot.polar = (random.random() * (self._stimulus_area_radius -
                                            dot.diameter / 2.0),
                         random.random() * TWO_PI)

            bad_position = False
            if not ignore_overlapping:
                # find bad_positions
                for c in self.dots:
                    if dot.distance(c) < self._min_gap:
                        bad_position = True
                        break  # for

            if not bad_position:
                return
            elif cnt > 3000:
                raise RuntimeError("Can't find a solution")

    def shuffle_all_positions(self, ignore_overlapping=False,
                              minimum_gap = None):
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        if minimum_gap is not None:
            self._min_gap = minimum_gap

        tmp_dots = self.dots
        self.dots = []
        for d in tmp_dots:
            self.randomize_dot_position(d, ignore_overlapping)
            self.dots.append(d)

    @property
    def minimum_gap(self):
        return self._min_gap

    @property
    def property_string(self):
        return self.property_names + "\n" + \
               str(self.properties).replace("[", "").replace("]", "")

    def get_array_csv_text(self, label=None, num_format="%10.2f",
                           variable_names=False):
        """Return the dot array as csv text

        Parameter
        ---------
        label : str or numeric, optional
            optional label of the dot array that is printed
        TODO
        variable_names : bool, optional
            if True variable name will be printed in the first line

        TODO add colour to csv
        """

        rtn = ""
        if variable_names:
            if label is not None:
                rtn += "label,"
            rtn += "cnt,x,y,diameter\n"

        for cnt, d in enumerate(self.dots):
            if label is not None:
                rtn += "{0},".format(label)
            rtn += "{0},".format(cnt + 1)
            rtn += num_format % d.x + "," + num_format % d.y + "," + \
                   num_format % d.diameter + "\n"
        return rtn

    def change_dot_colours(self, colour, dot_ids=None):
        """Change the colour of dots.

        Parameter
        ---------
        colour: colour
            the dot colour
        dot_ids: list of integers, optional
            if defined it set the colour of the dot defined by the id

        Note
        ----
        The function ignores dot limitations
        """

        if dot_ids is None:
            dot_ids = range(len(self.dots))
        for x in map(lambda i: self.dots[i], dot_ids):
            x.colour = colour

    def change_dot_pictures(self, picture, dot_ids=None):
        """Change the picture of dots.

        Parameter
        ---------
        colour: colour
            the dot colour
        dot_ids: list of integers, optional
            if defined it set the colour of the dot defined by the id

        """

        if dot_ids is None:
            dot_ids = range(len(self.dots))
        for x in map(lambda i: self.dots[i], dot_ids):
            x.picture = picture

    @property
    def property_names(self):
        return "n_dots,mean_dot_diameter,total_area,convex_hull_area,density, total_circumference"

    @property
    def properties(self):
        return [len(self.dots), self.mean_dot_diameter, self.total_area,
                self.convex_hull_area, self.density, self.total_circumference]

    @property
    def numpy_dot_positions(self):
        """list of tuples with xy coordinates"""
        return np.array(list(map(lambda x: x.xy, self.dots)))

    @property
    def numpy_dot_diameters(self):
        """list of tuples with xy coordinates"""
        return np.array(list(map(lambda x: x.diameter, self.dots)))


    @property
    def mean_dot_diameter(self):
        return self.numpy_dot_diameters.mean()

    @property
    def total_area(self):
        return sum(list(map(lambda x: x.area, self.dots)))

    @property
    def total_circumference(self):
        return sum(list(map(lambda x: x.circumference, self.dots)))

    @property
    def convex_hull_area(self):
        try:
            return ConvexHull(self.numpy_dot_positions).area
        except:
            return None

    @property
    def convex_hull_points(self):
        pos = self.numpy_dot_positions
        try:
            x = ConvexHull(pos).vertices
        except:
            return None

        return pos[x,:]

    @property
    def density_stimulus_area(self):
        """density takes into account the full possible stimulus area """
        return np.pi * self._stimulus_area_radius ** 2 / self.total_area

    @property
    def density(self):
        """density takes into account the convex hull"""
        try:
            return self.convex_hull_area / self.total_area
        except:
            return None

    @property
    def array_center(self):
        x = self.numpy_dot_positions
        minmax = np.array((np.min(x, axis=0), np.max(x, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    def center_dots(self):
        shift = self.array_center
        for d in self.dots:
            d.xy = (d.x - shift[0], d.y - shift[1])

    def fit_convex_hull_area(self, convex_hull_area, precision=0.01):
        """changes the convex hull area to a desired size with certain precision

        iterative method can takes some time.
        """

        current = self.convex_hull_area
        if current is None\
                or abs(convex_hull_area - current) < precision:
            return # It'S good already or not defined

        # iteratively determine scale
        scale = 1 # find good scale
        step = 0.1
        if convex_hull_area < current: # current too larger
            step *= -1

        # make a centered copy just to calc the radius
        centered_positions = NumpyPositionList(self.dots)
        centered_positions.xy -= self.array_center
        current = 0
        cnt = 0

        while abs(current - convex_hull_area) > precision:
            cnt += 1
            scale += step

            tmp = deepcopy(centered_positions)
            tmp.polar[:,0] *= scale
            tmp.calc_cartesian()
            current = ConvexHull(tmp.xy).area

            if (current < convex_hull_area and step < 0) or \
                    (current > convex_hull_area and step > 0):
                step *= -0.2 # change direction and finer grain

        # apply the scale
        for d in self.dots:
            d.pos_radius = d.pos_radius * scale

    def fit_density(self, density, precision=0.01, ratio_area_convex_hull_adaptation = 0.5):
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area and convex_hull (between 0 and 1)

        """

        # dens = convex_hull_area / total_area
        if ratio_area_convex_hull_adaptation<0 or ratio_area_convex_hull_adaptation>1:
            ratio_area_convex_hull_adaptation = 0.5

        density_change = density - self.density
        d_change_convex_hull = density_change * (1-ratio_area_convex_hull_adaptation)
        if abs(d_change_convex_hull) > 0:
            self.fit_convex_hull_area(convex_hull_area= self.total_area * (self.density+d_change_convex_hull),
                                  precision=precision)

        self.fit_total_area(total_area=self.convex_hull_area / density)



    def sort_dots_by_eccentricity(self):
        """Sort the dot list from inner to outer dots"""

        def dot_eccentricity(d):
            return d.pos_radius
        self.dots.sort(key=dot_eccentricity)

    def realign(self, minimum_gap = None):
        """Realigns the dots in order to remove all dots overlaps. If two dots
        overlap, the dots that is further apart from the arry center will be
        moved opposite to the direction of the other dot until there is no
        overlap (note: minimun_gap parameter). If two dots have exactly the same
        position the same position one is minimally shifted in a random direction.

        """

        if minimum_gap is not None:
            self._min_gap = minimum_gap


        d = NumpyPositionList(self.dots)
        d.jitter_identical_positions()
        shift_required = False

        # from inner to outer remove overlaps
        for i in np.argsort(d.polar[:, 0]):
            if d.remove_overlap_for_dot(dotid=i, min_gap=self._min_gap):
                shift_required = True

        # sqeeze in points that pop out of the stimulus area radius
        cnt = 0
        while True:
            too_far = np.where( (d.polar[:,0] + d.diameter//2) > self._stimulus_area_radius)[0] # find outlier
            if len(too_far)>0:
                # squeeze in outlier
                tmp = deepcopy(self.dots[too_far[0]])
                tmp.pos_radius = self._stimulus_area_radius - tmp.diameter//2 - 0.000000001
                d.xy[too_far[0], :] = tmp.xy
                # remove overlaps centered around new outlier position
                d.xy -= tmp.xy
                # remove all overlaps (inner to outer, i.e. starting with outlier)
                for i in np.argsort(d.polar[:, 0]):
                    d.remove_overlap_for_dot(dotid=i, min_gap=self._min_gap)

                # new pos for outlyer
                d.xy += tmp.xy # move back to old position
                shift_required = True
            else:
                break # end while loop

            cnt += 1
            if cnt>500:
                raise RuntimeError("Can't find solution when removing outlier (n=" + str(len(self.dots))+")")

        for c, xy in enumerate(d.xy):
            self.dots[c].xy = xy

        if not shift_required:
            return False
        else:
            self.realign() # recursion
            return True


    def fit_total_area(self, total_area):
        """dots will be realigned"""
        a_scale = (total_area / self.total_area)
        x = 2 / np.sqrt(np.pi)
        for d in self.dots:
            d.diameter = np.sqrt(d.area * a_scale) * x  # d=sqrt(4a/pi) = 2*sqrt(pi)*sqrt(a)

    def fit_mean_item_size(self, mean_dot_diameter):
        """dots will be realigned"""
        scale = mean_dot_diameter / self.mean_dot_diameter
        for d in self.dots:
            d.diameter = d.diameter * scale


    def number_deviant(self, change_numerosity):
        """number deviant

        Creates deviant with different numerosity but controlled for
        perceptual cues
        """

        # make a copy for the deviant
        deviant = self.copy()
        if len(self.dots)+change_numerosity <=0:
            deviant.dots = []
        else:
            # add or remove random dots
            for _ in range(abs(change_numerosity)):
                if change_numerosity<0:
                    # remove dots
                    deviant.dots.pop(random.randint(0, len(deviant.dots)-1))
                else:
                    # add
                    deviant.dots.append(deviant.create_random_dot(ignore_overlapping=False))

        return deviant