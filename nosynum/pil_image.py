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
                    background_colour=(255,255,255) ):
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


def write_pil_images_of_da_sequence(dot_array_sequence,
                      area_colour=None,
                      convex_hull_colour=None,
                      antialiasing=None,
                      background_colour=(255,255,255)):
    dot_array_sequence.images = []
    for da in dot_array_sequence.dot_arrays:
        im = create(dot_array=da, area_colour=area_colour,
                                      convex_hull_colour=convex_hull_colour,
                                      antialiasing=antialiasing,
                                      background_colour=background_colour)
        dot_array_sequence.images.append(im)


def pil_images_save(da_sequence_with_images, subfolder="stimuli", file_type="png",
                      file_prefix = "da-",
                      replace_images_by_filename = True):
    try:
        mkdir(subfolder)
    except:
        pass

    name = file_prefix + da_sequence_with_images.md5hash
    for x in range(len(da_sequence_with_images.images)):
        n = len(da_sequence_with_images.dot_arrays[x].dots)
        filename = unicode(path.join(subfolder, name + "-" + str(n) + "." + file_type))
        da_sequence_with_images.images[x].save(fp=filename, file_type=file_type)
        if replace_images_by_filename:
            da_sequence_with_images.images[x] = filename


class PILMakeDASequenceProcess(TemplateDASequenceProcess):

    def __init__(self, max_dot_array,
                 method=M_NO_FITTING,
                 n_trials=3,
                 auto_start_process=True,
                 save_images=False,
                 auto_delete_image_files = True,
                 subfolder="stimuli",
                 area_colour=(255, 255, 255),
                 convex_hull_colour=None,
                 antialiasing=None,
                 background_colour=(0,0,0)):
        """
        property: da_sequence, after processes finished
        Event(): sequence_available
        """

        super(PILMakeDASequenceProcess, self).__init__()

        self.max_dot_array = max_dot_array
        self.method = method
        self.n_trials = n_trials
        self.save_images = save_images
        self.subfolder = subfolder
        self.area_colour = area_colour
        self.convex_hull_colour = convex_hull_colour
        self.antialiasing = antialiasing
        self.background_colour = background_colour
        self.auto_delete_image_files = auto_delete_image_files

        if auto_start_process:
            self.start()

    def run(self):
        self._makeprocess = MakeDASequenceProcess(max_dot_array=self.max_dot_array,
                                                  method=self.method,
                                                  n_trials=self.n_trials,
                                                  auto_start_process=True)
        self._makeprocess.join()
        da_sequence = self._makeprocess.da_sequence
        write_pil_images_of_da_sequence( dot_array_sequence=da_sequence,
                          area_colour=self.area_colour,
                          convex_hull_colour=self.convex_hull_colour,
                          antialiasing=self.antialiasing,
                          background_colour=self.background_colour)
        if self.save_images:
            pil_images_save(da_sequence_with_images=da_sequence,
                        subfolder = self.subfolder,
                        replace_images_by_filename=True)
        da_sequence.auto_delete_image_files = self.auto_delete_image_files

        self.sequence_available.set()
        self._data_queue.put(da_sequence)
