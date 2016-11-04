"""
Dot Array
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math
import random
import pickle
from hashlib import sha1
import numpy as np
from nosynum import Dot, random_beta, shape_parameter_beta

def load_dot_array(filename):
    with open(filename, 'rb')as fl:
        rtn = pickle.load(fl)
    return rtn

class DotArray(object):
    def __init__(self, stimulus_area_radius, n_dots,
                 dot_diameter_mean,
                 dot_diameter_range=None,
                 dot_diameter_std=None, min_gap=1, dot_colour=None):

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

        self._min_gap = min_gap
        self._stimulus_area_radius = stimulus_area_radius
        self._dot_diameter_range = dot_diameter_range
        self._dot_diameter_mean = dot_diameter_mean
        # clear that this are expended values
        if dot_diameter_std <= 0:
            dot_diameter_std = None
        self._dot_diameter_std = dot_diameter_std
        self._dot_colour = dot_colour
        self._create_dots(n_dots)

    def _create_dots(self, n_dots):
        self._dots = []
        if self._dot_diameter_range is None or \
                        self._dot_diameter_std is None:
            # constant mean
            self._dots = [Dot(diameter=self._dot_diameter_mean,
                              colour=self._dot_colour) \
                          for _ in range(n_dots)]
        else:
            # draw diameter from beta distribution
            parameter = shape_parameter_beta(self._dot_diameter_range,
                                             self._dot_diameter_mean, self._dot_diameter_std)
            self._dots = [Dot(diameter=random_beta(
                self._dot_diameter_range, parameter),
                              colour=self._dot_colour) for _ in range(n_dots)]
        self._dot_limitation = None
        self.mix_positions()

    def mix_positions(self):
        # find new position for each dot
        # mixes always all position (ignores dot limitation)
        tmp_dots = self._dots
        self._dots = []
        for d in tmp_dots:
            cnt = 0
            while (True):
                cnt += 1
                polar = (random.random() * (self._stimulus_area_radius - \
                                            d.diameter / 2.0), random.random() * 2 * math.pi)
                d.xy = (polar[0] * math.cos(polar[1]),
                        polar[0] * math.sin(polar[1]))
                bad_position = False
                for c in self._dots:
                    if d.distance(c) < self._min_gap:
                        bad_position = True
                        break  # for
                if not (bad_position):
                    self._dots.append(d)
                    break  # while
                elif cnt > 10000:
                    raise RuntimeError("Can't find a solution")

    @property
    def hash_id(self):
        """secure hash (sha1) of the stimulus

        This is a unique id of this particular stimulus

        Notes
        -----
        Hash id is based on all available dots. That is, dot_limitations will
        be ignored.

        """

        return sha1(pickle.dumps(self._dots)).hexdigest()

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

        """  # TODO add colour to csv

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

    @property
    def dots(self):
        """list of Dots"""
        if self._dot_limitation is None:
            return self._dots
        else:
            return self._dots[:self._dot_limitation]

    def limit_number_of_dot(self, value):
        """Dot limitations

        If dot_limitations (integer) is defined (i.e. is not None), only the first
        x dots are considered in all subsequent methods and for all properties.
        Dot limitation will be reset with create_dots.  This function is, for
        instance, useful when displaying incrementally a dot pattern.

        if value is None dots array is unlimited.
        """

        try:
            x = int(value)
            if x < 0:
                x = 0
            elif x > len(self._dots):
                x = len(self._dots)
        except:
            x = None
        self._dot_limitation = x

    def set_dot_colour(self, colour, dot_ids=None):
        """Set the colour of dots.

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
            dot_ids = range(len(self._dots))
        for x in map(lambda x: self._dots[x], dot_ids):
            x.colour = colour

    @property
    def property_names(self):
        return "n_dots,mean_dot_diameter,total_area,convex_hull_area,density, total_circumference"

    @property
    def properties(self):
        return [len(self.dots), self.mean_dot_diameter, self.total_area,
                self.convex_hull_area, self.density, self.total_circumference]

    @property
    def mean_dot_diameter(self):
        return np.mean(map(lambda x: x.diameter, self.dots))

    @property
    def total_area(self):
        return sum(map(lambda x: x.area, self.dots))

    @property
    def total_circumference(self):
        return sum(map(lambda x: x.circumference, self.dots))

    @property
    def points(self):
        """list of tuples with xy coordinates"""
        return map(lambda x: x.xy, self.dots)

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

        points = sorted(self.points)
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
                self.mix_positions()
                remix += 1
            else:
                remix = 0
                if convex_hull_area > actual_cha:
                    self._stimulus_area_radius += 10
                else:
                    self._stimulus_area_radius -= 10
                self._create_dots(n_dots=len(self._dots))
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
                self.mix_positions()
                remix += 1
            else:
                remix = 0
                if density > actual_density:
                    self._stimulus_area_radius += 10
                else:
                    self._stimulus_area_radius -= 10
                self._create_dots(n_dots=len(self._dots))
            actual_density = self.density
            diff = abs(density - actual_density)

    def create_expyriment_stimulus(self, area_colour=None,
                                   convex_hull_colour=None, antialiasing=None):
        from expyriment.stimuli import Canvas, Circle, Line

        canvas = Canvas(size=[self._stimulus_area_radius * 2] * 2)
        if area_colour is not None:
            Dot(diameter=self._stimulus_area_radius,
                colour=area_colour).plot(canvas)
        if convex_hull_colour is not None:
            # plot convey hull
            hull = self.convex_hull
            hull.append(hull[0])
            last = None
            for p in hull:
                if last is not None:
                    Line(start_point=last, end_point=p, line_width=2,
                         colour=convex_hull_colour).plot(canvas)
                last = p
        map(lambda d: Circle(diameter=2 * int(d.radius), colour=d.colour,
                             line_width=0, position=d.xy).plot(canvas),
            self.dots)
        return canvas

    def create_pil_image(self, area_colour=None, convex_hull_colour=None,
                         antialiasing=True):
        """returns pil image"""
        import Image, ImageDraw

        pict_size = int(round(self._stimulus_area_radius * 2))

        def convert_pos(xy):
            j = int(float(pict_size) / 2)
            return (int(xy[0]) + j, int(-1 * xy[1]) + j)

        def draw_dot(draw, dot):
            if dot.colour is None:
                colour = (200, 200, 200)
            else:
                colour = dot.colour
            r = int(dot.radius)
            x, y = convert_pos(dot.xy)
            draw.ellipse((x - r, y - r, x + r, y + r), fill=colour)

        img = Image.new("RGB", (pict_size, pict_size), "black")
        draw = ImageDraw.Draw(img)
        if area_colour is not None:
            draw_dot(draw, Dot(x=0, y=0, diameter=self._stimulus_area_radius * 2),
                     colour=area_colour)
        if convex_hull_colour is not None:
            # plot convey hull
            hull = self.convex_hull
            hull.append(hull[0])
            last = None
            for p in hull:
                if last is not None:
                    draw.line(convert_pos(last) + convert_pos(p),
                              width=2, fill=convex_hull_colour)
                last = p
        map(lambda d: draw_dot(draw, d), self.dots)
        if antialiasing:
            img = img.resize((pict_size * 2, pict_size * 2))
            img = img.resize((pict_size, pict_size), Image.ANTIALIAS)
        return img

    def save_image(self, filename, file_type="PNG", area_colour=None,
                   convex_hull_colour=None, antialiasing=True):
        """Save"""
        img = self.create_pil_image(area_colour=area_colour,
                                    convex_hull_colour=convex_hull_colour,
                                    antialiasing=antialiasing)
        img.save(filename, file_type)

    def save_incremental_images(self, name, file_type="PNG", area_colour=None,
                                convex_hull_colour=None, antialiasing=True,
                                property_num_format="%10.2f"):
        """Saving incrementally.

        Each numerosity will be saved in a separate file by adding dots
        incrementally.

        returns
        -------
        rtn : string
            property list as string
        """

        old_dlim = self._dot_limitation
        properties_str = ""
        hash = self.hash_id
        for n_dots in range(1, len(self._dots) + 1):
            filename = name + "_incre_{0}.{1}".format(n_dots, file_type.lower())
            # print filename
            self.limit_number_of_dot(n_dots)
            self.save_image(filename=filename, file_type=file_type,
                            area_colour=area_colour,
                            convex_hull_colour=convex_hull_colour,
                            antialiasing=antialiasing)

            probs = self.properties
            properties_str += "{0},{1},{2},{3}".format(hash, filename, probs[0],
                                                       probs[1])
            for x in probs[2:]:
                if x is not None:
                    properties_str += "," + property_num_format % x
                else:
                    properties_str += ",None"
            properties_str += "\n"

        self.limit_number_of_dot(old_dlim)
        return "hash_id, filename, " + self.property_names + "\n" + properties_str
