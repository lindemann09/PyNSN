from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from os import path, mkdir
from multiprocessing import Pool
import codecs
from PIL import Image, ImageDraw
from . import Dot, DotArray
from .dot_array import PROPERTY_FILE_SUFFIX

def create_image(dot_array, area_colour=None,
                 convex_hull_colour=None,
                 antialiasing=True,
                 background_colour="white"):
    """returns pil image"""
    pict_size = int(round(dot_array._stimulus_area_radius * 2))
    def convert_pos(xy):
        j = int(float(pict_size) / 2)
        return (int(xy[0]) +j, int(- 1 *xy[1]) +j)

    def draw_dot(img, dot):
        if dot.colour is None:
            colour = (200 ,200 ,200)
        else:
            colour = dot.colour
        r = int(dot.radius)
        x ,y = convert_pos(dot.xy)
        if dot.picture is not None:
            pict = Image.open(dot.picture, "r")
            img.paste(pict, ( x -r , y -r))
        else:
            ImageDraw.Draw(img).ellipse(( x -r, y- r, x + r, y + r), fill=colour)

    img = Image.new("RGBA", (pict_size, pict_size), color=background_colour)
    if area_colour is not None:
        draw_dot(img, Dot(x=0, y=0, diameter=dot_array._stimulus_area_radius * 2,
                          colour=area_colour))
    if convex_hull_colour is not None:
        # plot convey hull
        hull = dot_array.convex_hull
        hull.append(hull[0])
        last = None
        draw = ImageDraw.Draw(img)
        for p in hull:
            if last is not None:
                draw.line(convert_pos(last) + convert_pos(p),
                          width=2, fill=convex_hull_colour)
            last = p
    list(map(lambda d: draw_dot(img, d), dot_array.dots))
    if antialiasing:
        img = img.resize((pict_size * 2, pict_size * 2), Image.ANTIALIAS)
        img = img.resize((pict_size, pict_size), Image.ANTIALIAS)
    return img

def save_incremental(dot_array,
                     name,
                     picture_folder=None,
                     area_colour=None,
                     convex_hull_colour=None, antialiasing=True,
                     property_num_format="%10.2f",
                     background_colour = "white",
                     save_properties =True):

    """Saving incrementally.
    Each numerosity will be saved in a separate file by adding dots
    incrementally.
    returns
    -------
    rtn : string
        property list as string
    """
    if picture_folder is not None:
        try:
            mkdir(picture_folder)
        except:
            pass
        name = path.join(picture_folder, name)

    parameter = map(lambda x: [dot_array, name, x, area_colour,
                               convex_hull_colour, antialiasing, property_num_format,
                               background_colour],
                    range(len(dot_array.dots) + 1))
    properties = Pool().map(_map_fnc_save_incremental, parameter)

    if save_properties:
        with codecs.open(name + PROPERTY_FILE_SUFFIX, "w", encoding="utf-8") as fl:
            fl.write("filename, " + dot_array.property_names + "\n")
            for p in properties:
                fl.write(p)


def create_and_save_multiple_incremental_dot_arrays(picture_folder="arrays",
                                                    area_radius=200,
                                                    area_colour=(200, 200, 200),
                                                    max_n_dots=20,
                                                    dot_diameter_range=None,
                                                    dot_diameter_std=None,
                                                    dot_colour=(0, 255, 0),
                                                    dot_size=20,
                                                    dot_picture=None,
                                                    n_dot_arrays=3,  # max 25
                                                    background_colour_pil="black"):

    """  TODO   """

    if n_dot_arrays > 24:
        raise RuntimeError("Max number of sets is 25")

    sets = list(map(lambda x: chr(x), range(ord("a"), ord("a") + n_dot_arrays)))

    txt = ""
    for s in sets:
        da = DotArray(stimulus_area_radius=area_radius,
                             n_dots=max_n_dots,
                             dot_diameter_mean=dot_size,
                             dot_diameter_range=dot_diameter_range,
                             dot_diameter_std=dot_diameter_std,
                             dot_colour=dot_colour,
                             dot_picture=dot_picture)

        save_incremental(dot_array=da,
                         area_colour=area_colour,
                         name="array_" + s,
                         picture_folder = picture_folder,
                         antialiasing=True,
                         background_colour = background_colour_pil)


def _map_fnc_save_incremental(parameter):
    # helper function for Pool().map()
    file_type = "PNG"
    (dot_array, name, n_dots, area_colour, convex_hull_colour,
                antialiasing, property_num_format, bkg_colour) = parameter
    filename = name + "_incre_{0}.{1}".format(n_dots, file_type.lower())
    print(filename)

    da = dot_array.copy(range(n_dots))
    pil_img = create_image(dot_array=da,
                           area_colour=area_colour,
                           convex_hull_colour=convex_hull_colour,
                           antialiasing=antialiasing,
                           background_colour=bkg_colour)
    pil_img.save(filename, file_type)

    property_txt = "{0},{1},{2}".format(filename,
                da.properties[0], da.properties[1])
    for x in da.properties[2:]:
        if x is not None:
            property_txt += "," + property_num_format%x
        else:
            property_txt += ",None"
    property_txt += "\n"
    return property_txt