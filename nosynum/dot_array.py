"""
Dot Array
"""
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import types
import math
import random
import pickle
from copy import deepcopy, copy
from hashlib import sha1
import numpy as np
from . import Dot, random_beta, shape_parameter_beta

TWO_PI = 2*math.pi

def load_dot_array(filename):
    with open(filename, 'rb')as fl:
        rtn = pickle.load(fl)
    return rtn

class DotArray(object):
    def __init__(self, stimulus_area_radius, n_dots,
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
            if isinstance(subset_dot_ids, types.IntType):
                subset_dot_ids = [subset_dot_ids]
            rtn.dots = []
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

        if minimum_gap is None:
            minimum_gap = self._min_gap

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
                    if dot.distance(c) < minimum_gap:
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
    def hash_id(self):
        """secure hash (sha1) of the stimulus

        This is a unique id of this stimulus

        Notes
        -----
        Hash id is based on all available dots.

        """

        return sha1(pickle.dumps(self.dots)).hexdigest()

    def __str__(self):
        return "dot pattern: \n" + self.property_names + "\n" + \
               str(self.properties).replace("[", "").replace("]", "")

    def get_array_csv_text(self, label=None, num_format="%10.2f",
                           variable_names=False, print_hash_id=False):
        """Return the dot array as csv text

        Parameter
        ---------
        label : str or numeric, optional
            optional label of the dot array that is printed
        TODO
        variable_names : bool, optional
            if True variable name will be printed in the first line
        print_hash_id : boolean, optional
            if True unique hash will be printed (default: False)

        TODO add colour to csv
        """

        rtn = ""
        if variable_names:
            if print_hash_id:
                rtn += "hash_id,"
            if label is not None:
                rtn += "label,"
            rtn += "cnt,x,y,diameter\n"

        hash_id = self.hash_id
        for cnt, d in enumerate(self.dots):
            if print_hash_id:
                rtn += "{0},".format(hash_id)
            if label is not None:
                rtn += "{0},".format(label)
            rtn += "{0},".format(cnt + 1)
            rtn += num_format % d.x + "," + num_format % d.y + "," + \
                   num_format % d.diameter + "\n"
        return rtn

    def save(self, filename):
        """Saving dot array object. Use load_array to load save array"""
        with open(filename, 'wb')as fl:
            pickle.dump(self, fl)

    def change_dot_colour(self, colour, dot_ids=None):
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
        for x in map(lambda x: self.dots[x], dot_ids):
            x.colour = colour

    def change_dot_picture(self, picture, dot_ids=None):
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
        for x in map(lambda x: self.dots[x], dot_ids):
            x.picture = picture

    @property
    def property_names(self):
        return "n_dots,mean_dot_diameter,total_area,convex_hull_area,density, total_circumference"

    @property
    def properties(self):
        return [len(self.dots), self.mean_dot_diameter, self.total_area,
                self.convex_hull_area, self.density, self.total_circumference]

    @property
    def mean_dot_diameter(self):
        return np.mean(list(map(lambda x: x.diameter, self.dots)))

    @property
    def total_area(self):
        return sum(list(map(lambda x: x.area, self.dots)))

    @property
    def total_circumference(self):
        return sum(list(map(lambda x: x.circumference, self.dots)))

    @property
    def all_positions(self):
        """list of tuples with xy coordinates"""
        return list(map(lambda x: x.xy, self.dots))

    @property
    def convex_hull_area(self):
        """area extended by a stimulus or convex hull area

        returns the area of the closed polygon self.convex_hull

        """

        p = self.convex_hull
        if len(p) < 2:
            return None
        segments = zip(p, p[1:] + [p[0]])
        return 0.5 * abs(sum(x0 * y1 - x1 * y0 \
                             for ((x0, y0), (x1, y1)) in segments))

    @property
    def density2(self):
        """density takes into account the full possible stimulus area """
        return math.pi * self._stimulus_area_radius ** 2 / self.total_area

    @property
    def density(self):
        """density takes into account the convex hull"""
        h = self.convex_hull_area
        a = self.total_area
        if h is not None and a > 0:
            return self.convex_hull_area / self.total_area
        else:
            return None

    @property
    def convex_hull(self):
        """Return list of points (tuples) of the convex hull in
        counter-clockwise order.

        Implements Andrew's monotone chain algorithm. O(n log n) complexity.
        Credit: http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
        """

        def cross(o, a, b):
            # 2D cross product of OA & OB vectors,
            # i.e. z-component of their 3D cross product.
            # Returns a positive value, if OAB makes a counter-clockwise turn,
            # negative for clockwise turn, and zero if the points are collinear.
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

        points = sorted(self.all_positions)
        if len(points) <= 1:
            return points
        # Build lower hull
        lower = []
        for p in points:
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
                lower.pop()
            lower.append(p)
        # Build upper hull
        upper = []
        for p in reversed(points):
            while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
                upper.pop()
            upper.append(p)
        # Concatenation of the lower and upper hulls gives the convex hull.
        return lower[:-1] + upper[:-1]


    def adapt_convex_hull_area(self, convex_hull_area, precision=0.01):
        """this function changes the area and remixes to get a desired
        convex_hull_area
        precision in percent between 1 < 0
        """
        remix = 0
        actual_cha = self.convex_hull_area
        diff = abs(convex_hull_area - actual_cha)
        abs_precision = convex_hull_area * precision
        while (diff > abs_precision):
            if diff < (abs_precision * 2) and remix <= 10:  # try with remixing
                self.shuffle_all_positions()
                remix += 1
            else:
                remix = 0
                if convex_hull_area > actual_cha:
                    self._stimulus_area_radius += 10
                else:
                    self._stimulus_area_radius -= 10
                self._create_dots(n_dots=len(self.dots))
            actual_cha = self.convex_hull_area
            diff = abs(convex_hull_area - actual_cha)


    def adapt_density(self, density, precision=0.01):
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0
        """
        remix = 0
        actual_density = self.density
        diff = abs(density - actual_density)
        abs_precision = density * precision
        while (diff > abs_precision):
            if diff < (abs_precision * 2) and remix <= 10:  # try with remixing
                self.shuffle_all_positions()
                remix += 1
            else:
                remix = 0
                if density > actual_density:
                    self._stimulus_area_radius += 10
                else:
                    self._stimulus_area_radius -= 10
                self._create_dots(n_dots=len(self.dots))
            actual_density = self.density
            diff = abs(density - actual_density)

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

        diff = Dot()
        self.sort_dots_by_eccentricity()
        for a in range(len(self.dots)):
            for b in range(a+1, len((self.dots))):
                if self.dots[a].xy == self.dots[b].xy:
                    # jitter
                    diff.polar = (0.1, random.random() * TWO_PI)
                    self.dots[a].xy = (self.dots[a].x + diff.x,
                                       self.dots[a].y + diff.y)
                    return self.realign()

                diff.xy = (self.dots[b].x - self.dots[a].x,
                           self.dots[b].y - self.dots[a].y)
                m_gap = self._min_gap + (self.dots[a].diameter +
                                         self.dots[b].diameter) // 2
                if diff.pos_radius < m_gap:
                    diff.pos_radius += m_gap
                    self.dots[b].xy = (self.dots[a].x + diff.x,
                                       self.dots[a].y + diff.y)
                    return self.realign()
        return True


    def pcc_deviant(self, num_dots, match="total_area"):
        """perceptual cue controlled (pcc) deviant

        Creates deviant with different numerosity but controlled for
        perceptual cues

        """

        deviant = deepcopy(self)
        diff = num_dots - len(self.dots)
        if diff<0:
            for x in range(abs(diff)):
                deviant.dots.pop(random.randint(0, len(deviant.dots)))
        else:
            pass # TODO add_dots

        if match == "total_area":
            scale = np.sqrt((self.total_area*4)/(num_dots*np.pi)) / \
                    np.sqrt((self.total_area*4) / (len(self.dots) * np.pi))
            for d in deviant.dots:
                d.diameter = scale * d.diameter
            return deviant