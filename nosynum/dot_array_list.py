__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import os
import math
from multiprocessing import Pool
import numpy as np
from dot_array import DotArray
from expyriment.stimuli import Picture

def get_list_of_incremental_dot_arrays(subfolder="arrays",
                                       area_radius=200,
                                       area_colour=(200, 200, 200),
                                       max_n_dots=20,
                                       dot_diameter_range=None,
                                       dot_diameter_std=None,
                                       dot_colour=(0, 255, 0),
                                       dot_size=20,
                                       pictformat="png",
                                       dot_picture=None,
                                       n_sets=3,  # max 25
                                       position=(0, 0),
                                       background_colour_pil="black",
                                       background_stimulus_expyriment=None):

    """ CALL THIS FUNCTION
    returns the array with all stimuli
    """

    if n_sets > 24:
        raise RuntimeError("Max number od sets is 25")

    sets = list(map(lambda x: chr(x), range(ord("a"), ord("a") + n_sets)))
    if not os.path.isdir(subfolder):
        os.mkdir(subfolder)
        txt = ""
        for s in sets:
            dot_array = DotArray(stimulus_area_radius=area_radius,
                                 n_dots=max_n_dots,
                                 dot_diameter_mean=dot_size,
                                 dot_diameter_range=dot_diameter_range,
                                 dot_diameter_std=dot_diameter_std,
                                 dot_colour=dot_colour,
                                 dot_picture=dot_picture,
                                 background_colour_pil=background_colour_pil,
                                 background_stimulus_expyriment=background_stimulus_expyriment)
            property_str = dot_array.save_incremental_images(
                area_colour=area_colour,
                name=subfolder + os.path.sep + "array_" + s,
                file_type=pictformat, antialiasing=True)
            txt += property_str

        # properties
        with file(subfolder + os.path.sep + "properties.csv", 'w') as fl:
            for x in txt:
                fl.write(x)

    arrays = []
    filename = subfolder + os.path.sep + "array_{0}_incre_{1}." + pictformat
    for s in sets:
        array_list = []
        for x in range(max_n_dots + 1):
            array_list.append(Picture(filename.format(s, x), position=position))
        arrays.append(array_list)

    return arrays


#### property matched dot arrays (under development)

# map functions
def _map_fnc_dot_array_list_save(parameter):
    parameter[0].save_image(filename=parameter[5], file_type=parameter[1],
                area_colour=parameter[2], convex_hull_colour=parameter[3],
                antialiasing=parameter[4])

def _map_fnc_adapt_density(parameter):
    parameter[0].adapt_density(density=parameter[1], precision=parameter[2])
    return parameter[0]
def _map_fnc_adapt_convex_hull_area(parameter):
    parameter[0].adapt_convex_hull_area(convex_hull_area=parameter[1], precision=parameter[2])
    return parameter[0]


class DotArrayListMatched(object):

    def __init__(self, number_list, subfolder=""):
        self.number_list = number_list
        self._dot_arrays = []
        self.subfolder = subfolder

    def create_method(self, method, dot_diameter,
            dot_diameter_std=None, stimulus_area_factor = 10,
            min_gap=1, dot_colour=None):
        """Creating dot array list with one of three methods. see Izard & Dehaene"""

        max_stimulus_area_radius =  math.sqrt( stimulus_area_factor * \
                                    max(self.number_list) * dot_diameter**2/4) + 1

        precision = 0.01 #TODO
        density = 8
        max_num = max(self.number_list)
        self._dot_arrays = []
        for num in self.number_list:
            if method==1:
                # constant convex hull, total area, density
                # variable: diameter
                stimulus_area_radius = max_stimulus_area_radius
                dot_diameter = 2*stimulus_area_radius / \
                            math.sqrt(num*stimulus_area_factor)
            elif method==2:
                # constant diameter, convex hull
                # variable: total area
                stimulus_area_radius = max_stimulus_area_radius
            elif method==3:
                # constant diameter, density
                # variable convex hull, total area
                stimulus_area_radius =  math.sqrt(density *num) * dot_diameter/2

            else:
                raise RuntimeError("Unknown method")

            if dot_diameter_std is None:
                dot_diameter_range = None
            else:
                dot_diameter_range = (dot_diameter-3*dot_diameter_std,
                                      dot_diameter+3*dot_diameter_std)

            dp = DotArray(n_dots = num,
                    stimulus_area_radius = stimulus_area_radius,
                    dot_diameter_mean = dot_diameter,
                    dot_diameter_range = dot_diameter_range,
                    dot_diameter_std = dot_diameter_std,
                    min_gap=min_gap,
                    dot_colour=dot_colour)
            #filename
            if len(self.subfolder)>2:
                dp.filename = "{0}/".format(self.subfolder)
            else:
                dp.filename = ""
            dp.filename += "num{0}_{1}.png".format(num, method)
            self._dot_arrays.append(dp)

        ## adaptations
        if method in [1, 3]:
            parameter = map(lambda x:[x, density, precision], self._dot_arrays)
            self._dot_arrays = Pool().map(_map_fnc_adapt_density, parameter)
        elif method==2:
            cha = max(map(lambda x:x.convex_hull_area, self._dot_arrays))
            precision = 0.01
            parameter = map(lambda x:[x, cha, precision], self._dot_arrays)
            self._dot_arrays = Pool().map(_map_fnc_adapt_convex_hull_area, parameter)

    def get_arrays_csv_text(self, num_format="%10.2f"):
        """Dot array as csv text

        Parameter
        ---------
        index : str or numeric, optional
            optinal unique index of the dot array that is printed before
            each line
        TODO
        """
        for cnt, ar in enumerate(self._dot_arrays):
            if cnt==0:
                rtn = ar.get_array_csv_text(cnt, num_format=num_format,
                                            variable_names=True)
            else:
                rtn += ar.get_array_csv_text(cnt, num_format=num_format,
                                            variable_names=False)
        return rtn

    @property
    def dot_arrays(self):
        """getter for dot arrays (list of all dot arrays)"""
        return self._dot_arrays

    @property
    def property_names(self):
        return self._dot_arrays[0].property_names

    @property
    def properties(self):
        prop = []
        for cnt, ar in enumerate(self._dot_arrays):
            prop.append([cnt] + ar.properties)
        return prop

    def extend(self, other_dot_array_list):
        self._dot_arrays.extend(other_dot_array_list._dot_arrays)

    def analysis(self):
        data = np.array(self.properties)
        def cor(label, var1,var2):
            r = np.corrcoef(var1, var2)[0, 1]
            print "{0}: R2={1}".format(label, r**2)

        cor("mean dot diameter", data[:,1], data[:, 2] )
        cor("total_area", data[:,1], data[:, 3])
        cor("convex_hull_area", data[:,1], data[:, 4])
        cor("density", data[:,1], data[:, 5])
        cor("total_circumference", data[:,1], data[:, 6])

    def save_arrays_csv_text(self, filename, num_format="%10.2f"):
        """save dot arrays as text files
        TODO
        """

        with file(filename, 'w') as fl:
            fl.write(self.get_arrays_csv_text(num_format=num_format))

    def save_properties(self, filename,  num_format="%10.2f"):
        """save properties to text file includes variable names

        TODO (do not use numpy save, first to cols should be int)

        """

        with file(filename, 'w') as fl:
            fl.write("id, " + self.property_names + "\n")
            for s in self.properties:
                fl.write("{0},{1}".format(s[0], s[1]))
                for x in s[2:]:
                    if x is not None:
                        fl.write("," + num_format%x)
                    else:
                        fl.write(",None")
                fl.write("\n")

    def save_images(self, name, file_type="PNG", area_colour=None,
                convex_hull_colour=None, antialiasing=True):
        """
        save a dot array list as pictures
        using multiprocessing
        name with path
        """

        #set filenames
        for cnt, d in enumerate(self._dot_arrays):
            d.filename = name + "_{0}_num{1}.{2}".format(cnt, d.n_dots,
                            file_type.lower())

        #list with parameter
        parameter = map(lambda x:[x, file_type, area_colour, convex_hull_colour,
                antialiasing, x.filename],
                self._dot_arrays)
        Pool().map(_map_fnc_dot_array_list_save, parameter)
