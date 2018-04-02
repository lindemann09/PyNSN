"""
Dot Array
"""
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import random
from copy import deepcopy, copy
from hashlib import md5
import numpy as np
from scipy.spatial import ConvexHull
from .dot import Dot

TWO_PI = 2*np.pi

def polar2cartesian(polar):
    """polar is an 2d-array representing polar coordinates (radius, angle)"""
    polar = np.array(polar)
    return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                     polar[:, 0] * np.sin(polar[:, 1])]).T

def cartesian2polar(xy, radii_only=False):
    xy = np.array(xy)
    radii = np.hypot(xy[:, 0], xy[:, 1])
    if radii_only:
        return radii
    else:
        return np.array([radii, np.arctan2(xy[:, 1], xy[:, 0])]).T

def short_md5_hash(unicode, hash_length=8):
    return md5(unicode.encode('utf-8')).hexdigest()[:hash_length]


class DotList(object):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    def __init__(self, xy_positions = (), diameters = ()):

        self.xy = np.array(xy_positions)
        self.diameters = np.array(diameters)

    def clear(self):
        self.xy = np.array([[]])
        self.diameters = np.array([])

    def polar(self):
        return cartesian2polar(self.xy)

    def distances(self, xy, diameter):
        """Distances toward a single point (xy, diameter) """
        if len(self.xy)==0:
            return np.array([])
        else:
             return np.hypot(self.xy[:,0] - xy[0], self.xy[:, 1] - xy[1]) - \
                            ((self.diameters + diameter) / 2.0)
    @property
    def rounded_xy(self):
        """rounded to integer"""
        return np.array(np.round(self.xy), dtype=np.int32)

    @property
    def rounded_diameters(self):
        """rounded to integer"""
        return np.array(np.round(self.diameters), dtype=np.int32)

    @property
    def center_of_outer_positions(self):
        minmax = np.array((np.min(self.xy, axis=0), np.max(self.xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def center_of_mass(self):
        weighted_sum = np.sum(self.xy *self.diameters[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self.diameters)

    @property
    def surface_areas(self):
        return np.pi * (self.diameters ** 2) / 4.0

    @property
    def circumferences(self):
        return np.pi * self.diameters

    @property
    def convex_hull_positions(self):
        x = ConvexHull(self.xy).vertices
        return self.xy[x,:]

    @property
    def convex_hull(self):
        """this convex_hull takes into account the dot diameter"""
        idx = ConvexHull(self.xy).vertices
        center = self.center_of_outer_positions
        polar_centered = cartesian2polar(self.xy[idx, :] - center)
        polar_centered[:,0] = polar_centered[:,0] + (self.diameters[idx]/2)
        xy = polar2cartesian(polar_centered) + center
        return xy

    def _jitter_identical_positions(self, jitter_size=0.1):
        """jitters points with identical position"""

        for idx, ref_dot in enumerate(self.xy):
            identical = np.where(np.all(np.equal(self.xy, ref_dot), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != idx:
                        self.xy[x, :] -= polar2cartesian([[jitter_size, random.random() * TWO_PI]])[0]

    def remove_overlap_for_dot(self, dot_id, minimum_gap):
        """remove overlap for one point"""


        dist = self.distances(self.xy[dot_id, :], self.diameters[dot_id])

        shift_required = False
        idx = np.where(dist < minimum_gap)[0].tolist() # overlapping dot ids
        if len(idx)>1:
            idx.remove(dot_id) # don't move yourself
            if np.sum(np.all(self.xy[idx, ] == self.xy[dot_id, :], axis=1)) > 0: # check if there is an identical position
                self._jitter_identical_positions()

            tmp_polar = cartesian2polar(self.xy[idx, :] - self.xy[dot_id, :])
            tmp_polar[:, 0] = 0.000000001 + minimum_gap - dist[idx]  # determine movement size
            xy = polar2cartesian(tmp_polar)
            self.xy[idx, :] = np.array([self.xy[idx, 0] + xy[:, 0], self.xy[idx, 1] + xy[:, 1]]).T
            shift_required = True

        return shift_required

    @property
    def prop_mean_dot_diameter(self):
        return self.diameters.mean()

    @property
    def prop_total_surface_area(self):
        return np.sum(self.surface_areas)

    @property
    def prop_total_circumference(self):
        return np.sum(self.circumferences)

    @property
    def prop_convex_hull_area_positions(self):
        return ConvexHull(self.xy).area

    @property
    def prop_convex_hull_area(self):
        return ConvexHull(self.convex_hull).area

    @property
    def prop_numerosity(self):
        return len(self.xy)

    @property
    def prop_density(self):
        """density takes into account the convex hull"""
        try:
            return self.prop_convex_hull_area / self.prop_total_surface_area # todo: positions conved hull
        except:
            return None




class DotArray(DotList):

    def __init__(self, max_array_radius, minimum_gap=1):
        """Dot array is restricted to a certain area and can generate random dots and
            be realigned """

        DotList.__init__(self, xy_positions=(), diameters=())
        self.max_array_radius = max_array_radius
        self.colours = np.array([[]])
        self.pictures = np.array([])
        self.minimum_gap = minimum_gap

    def clear(self):
        DotList.clear(self)
        self.colours = np.array([[]])
        self.pictures = np.array([])

    def append(self, xy, diameter, colour=None, picture=None):
        """append a single dot using numpy array"""

        colour = Dot.convert_colour(colour)
        if colour is None:
            colour = [None]*3

        if len(self.xy)==0:
            self.xy = np.array([xy])
            self.colours = np.array([colour])
        else:
            self.xy = np.append(self.xy, [xy], axis=0)
            self.colours = np.append(self.colours, [colour], axis=0)

        self.diameters = np.append(self.diameters, diameter)
        self.pictures = np.append(self.pictures, picture)

    def append_dot(self, dot):
        self.append(xy=[dot.x, dot.y], diameter=dot.diameter,
                    colour=dot.colour, picture=dot.picture)

    def get_dot(self, index):
        return Dot(x=self.xy[index, 0], y=self.xy[index, 1], diameter=self.diameters[index],
                   picture=self.pictures[index], colour=self.colours[index,:],)

    def delete_dot(self, index):
        self.xy = np.delete(self.xy, index, axis=0)
        self.diameters = np.delete(self.diameters, index)
        self.colours = np.delete(self.colours, index, axis=0)
        self.pictures = np.delete(self.pictures, index)

    def copy(self, subset_dot_ids=None):
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if subset_dot_ids is None:
            return deepcopy(self)
        else:
            rtn = copy(self)
            rtn.clear()
            if isinstance(subset_dot_ids, int):
                subset_dot_ids = [subset_dot_ids]
            for x in subset_dot_ids:
                rtn.append_dot(self.get_dot(x))
            return rtn

    def change_colours(self, colour, subset_dot_ids=None):
        """ """

        colour = Dot.convert_colour(colour)
        if isinstance(subset_dot_ids, int):
            subset_dot_ids = [subset_dot_ids]
        elif subset_dot_ids is None:
            subset_dot_ids = range(len(self.colours))

        self.colours[subset_dot_ids, :] = colour

    def change_colours_random_dots(self, colours, random_select_ratios=None):
        """change color of a part of the dots

        Tip: inhibit automatic logging and log manually, to get colous
        """

        if len(colours) != len(random_select_ratios):
            raise RuntimeError("Please specifiy the same amunt of ratios and colours (rgb1, rgb2,..).")

        # check random select ratios
        try:
            tmp = sum(random_select_ratios)
        except:
            tmp = 0
        if tmp !=1:
            raise RuntimeError("Error: sum(random_select_ratios) has to be 1!")

        n_dots = np.array(np.round(self.prop_numerosity * np.array(random_select_ratios)), dtype=np.int)
        idx = list(range(self.prop_numerosity))
        random.shuffle(idx)
        for c, n in zip(colours[:-1], n_dots[:-1]):
            i = idx[:n]
            idx[:n] = []
            print(c)
            self.change_colours(colour=c, subset_dot_ids=i)

        self.change_colours(colour=colours[-1], subset_dot_ids=idx)



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
        for i in np.argsort(cartesian2polar(self.xy, radii_only=True)):
            if self.remove_overlap_for_dot(dot_id=i, minimum_gap=self.minimum_gap):
                shift_required = True

        # sqeeze in points that pop out of the stimulus area radius
        cnt = 0
        while True:
            radii= cartesian2polar(self.xy, radii_only=True)
            too_far = np.where((radii + self.diameters // 2) > self.max_array_radius)[0] # find outlier
            if len(too_far)>0:

                # squeeze in outlier
                polar = cartesian2polar([self.xy[too_far[0],:]])[0]
                polar[0] = self.max_array_radius - self.diameters[too_far[0]] // 2 - 0.000000001 # new radius
                new_xy = polar2cartesian([polar])[0]
                self.xy[too_far[0], :] = new_xy

                # remove overlaps centered around new outlier position
                self.xy -= new_xy
                # remove all overlaps (inner to outer, i.e. starting with outlier)
                for i in np.argsort(cartesian2polar(self.xy, radii_only=True)):
                    self.remove_overlap_for_dot(dot_id=i, minimum_gap=self.minimum_gap)
                # new pos for outlyer
                self.xy += new_xy # move back to old position
                shift_required = True
            else:
                break # end while loop

            cnt += 1
            if cnt>20:
                error = True
                break

        if error:
            return False, "Can't find solution when removing outlier (n=" + \
                                     str(self.prop_numerosity) + ")"
        if not shift_required:
            return True, ""
        else:
            return self.realign() # recursion

    @property
    def prop_density_max_field(self):
        """density takes into account the full possible dot area """
        return np.pi * self.max_array_radius ** 2 / self.prop_total_surface_area

    @property
    def properties(self):
        return [self.prop_numerosity, self.prop_mean_dot_diameter, self.prop_total_surface_area,
                self.prop_convex_hull_area, self.prop_density, self.prop_total_circumference]

    def get_property_string(self, variable_names=False):
        rtn = ""
        if variable_names:
            rtn += "object_id, " + ", ".join(self.property_names) + "\n"

        rtn += self.object_id + ","
        rtn += str(self.properties).replace("[", "").replace("]", "") + "\n"
        return rtn

    @property
    def property_names(self):
        return ("n_dots", "mean_dot_diameter", "total_surface_area", "convex_hull_area",
                "density", "total_circumference")

    @property
    def object_id(self):
        """md5_hash of csv (counter, position, diameter only)"""

        csv = self.get_csv(num_format="%7.2f", object_id_column=False,
                           variable_names=False, num_idx_column=False,
                           colour_column=False, picture_column=False)
        return short_md5_hash(csv)

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
                rtn += "object_id,"
            if num_idx_column:
                rtn += "num_id,"
            rtn += "x,y,diameter"
            if colour_column:
                rtn += ",r,g,b"
            if picture_column:
                rtn += ",file"
            rtn += "\n"

        if object_id_column:
            obj_id = self.object_id

        for cnt in range(len(self.xy)):
            if object_id_column:
                rtn += "{0},".format(obj_id)
            if num_idx_column:
                rtn += "{},".format(self.prop_numerosity)
            rtn += num_format % self.xy[cnt, 0] + "," + num_format % self.xy[cnt, 1] + "," + \
                   num_format % self.diameters[cnt]
            if colour_column:
                colour = self.colours[cnt]
                if colour is not None and len(colour)==3:
                    rtn += ",{},{},{}".format(colour[0], colour[1], colour[2])
                else:
                    rtn += ", None, None, None"
            if picture_column:
                rtn += ", {}".format(self.pictures[cnt])
            rtn += "\n"
        return rtn

    def get_all_dots(self):
        rtn = []
        for xy, dia, col, pict in zip(self.xy, self.diameters,
                                      self.colours, self.pictures):

            rtn.append(Dot(x=xy[0], y=xy[1], diameter=dia,
                           colour=col, picture=pict))
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
            proposal_xy = polar2cartesian([proposal_polar])[0]

            bad_position = False
            if not ignore_overlapping:
                # find bad_positions
                dist = self.distances(proposal_xy, dot_diameter)
                idx = np.where(dist < self.minimum_gap)[0] # overlapping dot ids
                bad_position = len(idx)>0

            if not bad_position:
                return proposal_xy
            elif cnt > 3000:
                raise RuntimeError("Can't find a solution")

    def shuffle_all_positions(self, ignore_overlapping=False): #fixme centering?
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        diameters = self.diameters
        self.diameters = np.array([])
        self.xy = np.array([])
        for d in diameters:
            xy = self.random_free_dot_position(d, ignore_overlapping=ignore_overlapping)
            self.diameters = np.append(self.diameters, d)
            if len(self.xy) == 0:
                self.xy = np.array([xy])
            else:
                self.xy = np.append(self.xy, [xy], axis=0)

    def number_deviant(self, change_numerosity):
        """number deviant
        """

        # make a copy for the deviant
        deviant = self.copy()
        if self.prop_numerosity + change_numerosity <=0:
            deviant.clear()
        else:
            # add or remove random dots
            for _ in range(abs(change_numerosity)):
                rnd = random.randint(0, deviant.prop_numerosity-1) # strange: do not use np.rand, because parallel running process produce the same numbers
                if change_numerosity<0:
                    # remove dots
                    deviant.delete_dot(rnd)
                else:
                    # add
                    d = self.get_dot(rnd) # copy a random dot
                    d.xy = self.random_free_dot_position(dot_diameter=d.diameter)
                    deviant.append_dot(d)

        return deviant

    def match_total_surface_area(self, surface_area):
        # changes diameter
        a_scale = (surface_area / self.prop_total_surface_area)
        self.diameters = np.sqrt(self.surface_areas * a_scale) * 2/np.sqrt(np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)


    def match_mean_dot_diameter(self, mean_dot_diameter):
        # changes diameter

        scale = mean_dot_diameter / self.prop_mean_dot_diameter
        self.diameters = self.diameters * scale

    def match_total_circumference(self, total_circumference):
        # linear to fit_mean_dot_diameter, but depends on numerosity
        mean_dot_diameter = total_circumference/ (self.prop_numerosity * np.pi)
        self.match_mean_dot_diameter(mean_dot_diameter)


    def match_convex_hull_area(self, convex_hull_area, precision=0.0001,
                               center_array = True):
        """changes the convex hull area to a desired size with certain precision

        iterative method can takes some time.
        """

        current = self.prop_convex_hull_area

        if current is None:
            return # not defined

        # iteratively determine scale
        scale = 1 # find good scale
        step = 0.1
        if convex_hull_area < current: # current too larger
            step *= -1

        # centered  points
        old_center = self.center_of_outer_positions
        centered_polar = cartesian2polar(self.xy - old_center)

        while abs(current - convex_hull_area) > precision:
            scale += step

            self.xy = polar2cartesian(centered_polar * [scale, 1])
            current = self.prop_convex_hull_area

            if (current < convex_hull_area and step < 0) or \
                    (current > convex_hull_area and step > 0):
                step *= -0.2 # change direction and finer grain

        if not center_array:
            self.xy += old_center

    def match_density(self, density, precision=0.01, ratio_convex_hull2area_adaptation = 0.5):
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """

        # dens = convex_hull_area / total_surface_area
        if ratio_convex_hull2area_adaptation<0 or ratio_convex_hull2area_adaptation>1:
            ratio_convex_hull2area_adaptation = 0.5

        area_change100 = (self.prop_convex_hull_area / density) - self.prop_total_surface_area
        d_change_area = area_change100 * (1 - ratio_convex_hull2area_adaptation)
        if abs(d_change_area) > 0:
            self.match_total_surface_area(surface_area=self.prop_total_surface_area + d_change_area)

        self.match_convex_hull_area(convex_hull_area=self.prop_total_surface_area * density,
                                        precision=precision)

