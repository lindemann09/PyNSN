from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from os import path, mkdir
from PIL import Image, ImageDraw

from . import Dot
from .dot_array_sequences import M_NO_FITTING
from .multi_processing import MakeDASequenceProcess, TemplateDASequenceProcess

def create(dot_array,
                    area_colour=None,
                    convex_hull_colour=None,
                    antialiasing=None,  #TODO
                    background_colour="#ffffff"):
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image"""

    pict_size = int(round(dot_array.definition.stimulus_area_radius * 2))
    def convert_pos(xy):
        j = int(float(pict_size) / 2)
        return (int(xy[0]) +j, int(- 1 * xy[1]) + j)


    def draw_dot(img, dot):
        if dot.colour is None:
            colour = (200, 200, 200)
        else:
            colour = dot.colour
        r = dot.diameter // 2
        x, y = convert_pos(dot.xy)
        if dot.picture is not None:
            pict = Image.open(dot.picture, "r")
            img.paste(pict, (x - r, y - r))
        else:
            ImageDraw.Draw(img).ellipse((x - r, y - r, x + r, y + r), fill=colour)


    img = Image.new("RGBA", (pict_size, pict_size), color=background_colour)
    if area_colour is not None:
        draw_dot(img, Dot(x=0, y=0, diameter=dot_array.definition.stimulus_area_radius * 2,
                          colour=area_colour))
    if convex_hull_colour is not None:
        # plot convey hull
        hull = dot_array.convex_hull_points
        hull = list(hull) + [hull[0]]
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


def create_and_save(dot_array, filename,
           file_type="PNG",
           area_colour=None,
           convex_hull_colour=None,
           antialiasing=None,
           background_colour="#ffffff"):
    """see create"""

    pil_img = create(dot_array=dot_array, area_colour=area_colour,
           convex_hull_colour=convex_hull_colour, antialiasing=antialiasing,
           background_colour=background_colour)
    pil_img.save(filename, file_type)


def dasequence2images(dot_array_sequence,
                      image_filename,
                      subfolder="stimuli",
                      image_file_type="png",
                      area_colour=None,
                      convex_hull_colour=None,
                      antialiasing=None,
                      background_colour="#ffffff"):
    try:
        mkdir(subfolder)
    except:
        pass

    tmp = path.splitext(image_filename)
    for da in dot_array_sequence.dot_arrays:
        im = create(dot_array=da, area_colour=area_colour,
                                      convex_hull_colour=convex_hull_colour,
                                      antialiasing=antialiasing,
                                      background_colour=background_colour)
        filename = path.join(subfolder, tmp[0] + "-" + str(len(da.dots)) + tmp[1] + "." + image_file_type)
        im.save(fp=filename, file_type=image_file_type)



class PILMakeDASequenceProcess(TemplateDASequenceProcess):

    def __init__(self, max_dot_array,
                 method=M_NO_FITTING,
                 n_trials=3,
                 auto_start_process=True,
                 image_filename="array",
                 subfolder="stimuli",
                 image_file_type="png",
                 area_colour=None,
                 convex_hull_colour=None,
                 antialiasing=None,
                 background_colour="#ffffff"):

        super(PILMakeDASequenceProcess, self).__init__()

        self.max_dot_array = max_dot_array
        self.method = method
        self.n_trials = n_trials
        self.image_filename = image_filename
        self.subfolder = subfolder
        self.image_file_type = image_file_type
        self.area_colour = area_colour
        self.convex_hull_colour = convex_hull_colour
        self.antialiasing = antialiasing
        self.background_colour = background_colour

        if auto_start_process:
            self.start()

    def run(self):
        self._makeprocess = MakeDASequenceProcess(max_dot_array=self.max_dot_array,
                                                  method=self.method,
                                                  n_trials=self.n_trials,
                                                  auto_start_process=True)
        self._makeprocess.join()
        dasequence2images( dot_array_sequence=self._makeprocess.da_sequence,
                          image_filename=self.image_file_type,
                          subfolder=self.subfolder,
                          image_file_type=self.image_file_type,
                          area_colour=self.area_colour,
                          convex_hull_colour=self.convex_hull_colour,
                          antialiasing=self.antialiasing,
                          background_colour=self.background_colour)

        self.sequence_available.set()
        self._data_queue.put(self._makeprocess.da_sequence)
