"""
Dot Array
"""
from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from os import listdir
import types
import math
from copy import deepcopy, copy
from hashlib import md5
import numpy as np
from scipy.spatial import ConvexHull
from . import random_beta
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

def my_md5_hash(unicode, hash_length=8):
    return md5(unicode.encode('utf-8')).hexdigest()[:hash_length]



class DotList(object):
    """Numpy Position list for optimized for numpy calculations


    Position + diameter
    """

    def __init__(self, xy_positions = (), diameters = ()):

        self.xy = np.array(xy_positions)
        self.diameters = np.array(diameters)

    def clear(self):
        self.xy = np.array([])
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
    def xy_rounded_to_int(self):
        return np.array(np.round(self.xy), dtype=np.int32)

    @property
    def diameter_rounded_to_int(self):
        return np.array(np.round(self.diameters), dtype=np.int32)

    @property
    def center_of_positions(self):
        minmax = np.array((np.min(self.xy, axis=0), np.max(self.xy, axis=0)))
        return np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)

    @property
    def convex_hull_points(self):
        x = ConvexHull(self.xy).vertices
        return self.xy[x,:]

    def remove_overlap_for_dot(self, dot_id, minimum_gap=0):
        """remove overlap for one point"""

        dist = self.distances(self.xy[dot_id, :], self.diameters[dot_id])

        shift_required = False
        idx = np.where(dist < minimum_gap)[0].tolist() # overlapping dot ids
        if len(idx)>1:
            idx.remove(dot_id) # don't move yourself
            tmp_polar = cartesian2polar(self.xy[idx, :] - self.xy[dot_id, :])
            tmp_polar[:, 0] = 0.000000001 + minimum_gap - dist[idx]  # determine movement size
            xy = polar2cartesian(tmp_polar)
            self.xy[idx, :] = np.array([self.xy[idx, 0] + xy[:, 0], self.xy[idx, 1] + xy[:, 1]]).T
            shift_required = True

        return shift_required

    @property
    def prop_mean_dot_diameter(self):
        return self.xy.mean()

    @property
    def prop_total_area(self):
        return np.sum(np.pi * (self.diameters**2) / 4.0)

    @property
    def prop_total_circumference(self):
        return np.sum(np.pi * self.diameters)

    @property
    def prop_convex_hull_area(self):
        return ConvexHull(self.xy).area

    @property
    def prop_numerosity(self):
        return len(self.xy)

    @property
    def prop_density(self):
        """density takes into account the convex hull"""
        try:
            return self.prop_convex_hull_area / self.prop_total_area
        except:
            return None

    #todo def sort_dots_by_eccentricity(self):
    #    """Sort the dot list from inner to outer dots"""
    #
    #    def dot_eccentricity(d):
    #        return d.pos_radius
    #    self.dots.sort(key=dot_eccentricity)

class NumpyDotArray(DotList):

    def __init__(self, max_array_radius, xy_positions=(), diameters=(), colours = (),
                            pictures=(), minimum_gap=1):
        """Dot array is restricted to a certain area and can generate random dots and
            be realigned """

        DotList.__init__(self, xy_positions=xy_positions, diameters=diameters)
        self.max_array_radius = max_array_radius
        self.colours = np.array(colours)
        self.pictures = np.array(pictures)

    def clear(self):
        DotList.clear(self)
        self.colours = np.array([])
        self.pictures = np.array([])

    def append(self, xy, diameter, colour=None, picture=None):
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
        np.delete(self.xy, index, axis=0)
        np.delete(self.diameters, index)
        np.delete(self.colours, index, axis=0)
        np.delete(self.pictures, index)

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

    def _jitter_identical_positions(self, jitter_size=0.1):
        """jitters points with identical position"""

        for idx, ref_dot in enumerate(self.xy):
            identical = np.where(np.all(np.equal(self.xy, ref_dot), axis=1))[0]  # find identical positions
            if len(identical) > 1:
                for x in identical:  # jitter all identical positions
                    if x != idx:
                        self.xy[x, :] -= polar2cartesian([[jitter_size, np.random.random() * TWO_PI]])[0]

    def realign(self, minimum_gap):
        """Realigns the dots in order to remove all dots overlaps. If two dots
        overlap, the dots that is further apart from the arry center will be
        moved opposite to the direction of the other dot until there is no
        overlap (note: minimun_gap parameter). If two dots have exactly the same
        position the same position one is minimally shifted in a random direction.

        """

        self._jitter_identical_positions()
        shift_required = False
        outlier_error = False

        # old_xy = np.copy(self.xy) # keep a copy in the case no solution is found

        # from inner to outer remove overlaps
        for i in np.argsort(cartesian2polar(self.xy, radii_only=True)):
            if self.remove_overlap_for_dot(dot_id=i, minimum_gap=minimum_gap):
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
                    self.remove_overlap_for_dot(dot_id=i, minimum_gap=minimum_gap)
                # new pos for outlyer
                self.xy += new_xy # move back to old position
                shift_required = True
            else:
                break # end while loop

            cnt += 1
            if cnt>10:
                outlier_error = True
                break

        if outlier_error:
            raise RuntimeError("Can't find solution when removing outlier (n=" + \
                                    str(self.prop_numerosity) + ")")
        if not shift_required:
            return False
        else:
            self.realign(minimum_gap=minimum_gap) # recursion
            return True

    @property
    def prop_density_max_array(self):
        """density takes into account the full possible dot area """
        return np.pi * self.max_array_radius ** 2 / self.prop_total_area

    @property
    def properties(self):
        return [self.prop_numerosity, self.prop_mean_dot_diameter, self.prop_total_area,
                self.prop_convex_hull_area, self.prop_density, self.prop_total_circumference]

    def get_property_string(self, varnames=False):
        rtn = ""
        if varnames:
            rtn += "hash, " + ", ".join(self.property_names) + "\n"

        rtn += self.md5hash + ","
        rtn += str(self.properties).replace("[", "").replace("]", "")
        return rtn

    @property
    def property_names(self):
        return ("n_dots", "mean_dot_diameter", "total_area", "convex_hull_area",
                "density", "total_circumference")

    @property
    def md5hash(self):
        """md5_hash of csv (counter, position, diameter only)"""

        csv = self.get_csv(num_format="%7.2f", hash_column=False,
                           variable_names=False, n_dots_column=False,
                           colour_column=False, picture_column=False)
        return my_md5_hash(csv)

    def get_csv(self, num_format="%7.2f", hash_column=False,
                variable_names=True, n_dots_column=False,
                colour_column=False, picture_column=False):
        """Return the dot array as csv text

        Parameter
        ---------
        variable_names : bool, optional
            if True variable name will be printed in the first line

        """

        rtn = ""
        if variable_names:
            if hash_column:
                rtn += "hash,"
            if n_dots_column:
                rtn += "n_dots,"
            rtn += "dot_cnt,x,y,diameter"
            if colour_column:
                rtn += ",r,g,b"
            if picture_column:
                rtn += ",file"
            rtn += "\n"

        if hash_column:
            hash = self.md5hash

        for cnt in range(len(self.xy)):
            if hash_column:
                rtn += "{0},".format(hash)
            if n_dots_column:
                rtn += "{},".format(self.prop_numerosity)
            rtn += "{0},".format(cnt + 1)
            rtn += num_format % self.xy[cnt, 0] + "," + num_format % self.xy[cnt, 0] + "," + \
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


    @property
    def dots(self):  # fixme: to be compatible with old dot array
        rtn = []
        for xy, dia, col, pict in zip(self.xy, self.diameters,
                                      self.colours, self.pictures):

            rtn.append(Dot(x=xy[0], y=xy[1], diameter=dia,
                           colour=col, picture=pict))
        return rtn



class NumpyRandomDotArray(NumpyDotArray):

    def __init__(self, n_dots,
                 dot_array_definition):

        """Create a Random Dot Array
        """
        self.definition = copy(dot_array_definition)
        NumpyDotArray.__init__(self, max_array_radius= self.definition.stimulus_area_radius)
        self._create_dots(n_dots=n_dots)


    def _create_dots(self, n_dots):
        for _ in range(n_dots):
            xy, dia = self.create_random_dot(ignore_overlapping=False)
            self.append(xy=xy, diameter=dia,
                        colour=self.definition.dot_colour,
                        picture=self.definition.dot_picture)

    def create_random_dot(self, ignore_overlapping=False):
        """returns a random (xy, diameter)

        if no available space can be found it raise an Runtime Error
        """

        if self.definition.dot_diameter_range is None or \
                        self.definition.dot_diameter_std is None:
            # constant mean
            diameter=self.definition.dot_diameter_mean
        else:
            # draw diameter from beta distribution
            parameter = random_beta.shape_parameter_beta(self.definition.dot_diameter_range,
                                             self.definition.dot_diameter_mean,
                                             self.definition.dot_diameter_std)
            diameter=random_beta.random_beta(
                self.definition.dot_diameter_range, parameter)

        xy = self.get_random_dot_position(dot_diameter=diameter,
                                         ignore_overlapping=ignore_overlapping)

        return xy, diameter

    def get_random_dot_position(self, dot_diameter,
                               ignore_overlapping=False,
                               minimum_gap=None):
        """moves a dot to an available random position

        raise exception if not found"""

        if minimum_gap is None:
            minimum_gap = self.definition.minium_gap

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
                idx = np.where(dist < minimum_gap)[0] # overlapping dot ids
                bad_position = len(idx)>0

            if not bad_position:
                return proposal_xy
            elif cnt > 3000:
                raise RuntimeError("Can't find a solution")

    def shuffle_all_positions(self, ignore_overlapping=False,
                              minimum_gap = None):
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        diameters = self.diameters
        self.diameters = np.array([])
        self.xy = np.array([])
        for d in diameters:
            xy = self.get_random_dot_position(d, ignore_overlapping=ignore_overlapping,
                                        minimum_gap=minimum_gap)
            self.diameters = np.append(self.diameters, d)
            if len(self.xy) == 0:
                self.xy = np.array([xy])
            else:
                self.xy = np.append(self.xy, [xy], axis=0)

