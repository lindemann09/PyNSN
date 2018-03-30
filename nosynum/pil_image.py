from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from os import path, mkdir
from PIL import Image, ImageDraw
import numpy as np

from .dot_array_sequences import DASequence, M_NO_FITTING
from .multi_processing import TemplateDASequenceProcess


def _draw_dot(img, xy, diameter, colour=None, picture=None):
    # draw a dot on an image

    if colour is None:
        colour = (255, 255, 255)
    else:
        colour = colour

    r = diameter // 2
    if picture is not None:
        pict = Image.open(picture, "r")
        img.paste(pict, (xy[0]-r, xy[1]-r))
    else:
        ImageDraw.Draw(img).ellipse((xy[0]-r, xy[1]-r, xy[0]+r, xy[1]+r), fill=tuple(colour))

def _convert_pos(xy, image_size):
    """convert dot pos to pil image coordinantes"""
    return (xy * [1,-1]) + image_size//2

def _draw_convex_hull(img, convex_hull, convex_hull_colour):
    # plot convey hull
    hull = np.append(convex_hull, [convex_hull[0]], axis=0)
    last = None
    draw = ImageDraw.Draw(img)
    for p in hull:
        if last is not None:
            draw.line(np.append(last, p).tolist(),
                      width=2, fill=convex_hull_colour)
        last = p


def create(dot_array,
           colour_area=None,
           colour_convex_hull_positions=None,
           colour_convex_hull_dots=None,
           colour_center_of_mass = None,
           colour_center_of_outer_positions=None,
           antialiasing=None,  #TODO
           colour_background=(0, 0, 0)):
    """use PIL colours (see PIL.ImageColor.colormap)

    returns pil image"""

    image_size = int(round(dot_array.max_array_radius * 2))
    img = Image.new("RGBA", (image_size, image_size), color=colour_background)
    if colour_area is not None:
        _draw_dot(img, xy=_convert_pos(np.zeros(2), image_size),
                  diameter=image_size, colour=colour_area)

    for xy, d, c in zip(_convert_pos(dot_array.rounded_xy, image_size),
                        dot_array.rounded_diameters,
                        dot_array.colours):
        _draw_dot(img, xy=xy, diameter=d, colour=c)

    if colour_convex_hull_positions is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull = _convert_pos(dot_array.convex_hull_positions, image_size),
                          convex_hull_colour=colour_convex_hull_positions)

    if colour_convex_hull_dots is not None:
        # plot convey hull
        _draw_convex_hull(img=img,
                          convex_hull = _convert_pos(dot_array.convex_hull_dots, image_size),
                          convex_hull_colour=colour_convex_hull_dots)

    if colour_center_of_mass is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_mass, image_size),
                  diameter=10, colour=colour_center_of_mass)
    if colour_center_of_outer_positions is not None:
        _draw_dot(img, xy=_convert_pos(dot_array.center_of_outer_positions, image_size),
                  diameter=10, colour=(0, 255,0))

    if antialiasing:
        img = img.resize((image_size * 2, image_size * 2), Image.ANTIALIAS)
        img = img.resize((image_size, image_size), Image.ANTIALIAS)
    return img


def write_pil_images_of_da_sequence(dot_array_sequence,
                      area_colour=None,
                      convex_hull_colour=None,
                      antialiasing=None,
                      background_colour=(255,255,255)):
    """note: rounds dot array to intergers """
    dot_array_sequence.images = []
    for da in dot_array_sequence.dot_arrays:
        da.round_dot_paramter_to_integer()
        im = create(dot_array=da, colour_area=area_colour,
                    convex_hull_colour=convex_hull_colour, #fixme colours convex_hulls and center
                    antialiasing=antialiasing,
                    colour_background=background_colour)
        dot_array_sequence.images.append(im)


def pil_images_save(da_sequence_with_images, folder, file_type="png",
                    file_prefix = "da-",
                    replace_images_by_filename = True):
    try:
        mkdir(folder)
    except:
        pass

    name = file_prefix + da_sequence_with_images.md5hash
    for x in range(len(da_sequence_with_images.images)):
        n = len(da_sequence_with_images.dot_arrays[x].dots)
        filename = path.join(folder, name + u"-" + str(n) + u"." + file_type)
        da_sequence_with_images.images[x].save(fp=filename, file_type=file_type)
        if replace_images_by_filename:
            da_sequence_with_images.images[x] = filename


class PILMakeDASequenceProcess(TemplateDASequenceProcess):
    # making images in seperate prcess
    # DO NOT USE SAVE_IMAGES=FALSE, BECAUSE OF LARGE MEMORY REQUIRMENTS
    # AND TRANSFER OF PIL IMAGES BETWEEN PROCESSES IS SLOW
    # use multiprocessing.MakeDASequenceProcess and write_pil_images_of_da_sequence
    # in main process.

    def __init__(self, max_dot_array,
                 min_numerosity,
                 method=M_NO_FITTING,
                 n_trials=3,
                 auto_delete_image_files = True,
                 folder="tmp_stimuli",
                 area_colour=(255, 255, 255),
                 convex_hull_colour=None,
                 antialiasing=None,
                 background_colour=(0,0,0),
                 sqeeze_factor=None,
                 save_images=True):
        """
        property: da_sequence, after processes finished
        Event(): sequence_available
        """

        super(PILMakeDASequenceProcess, self).__init__()

        self.max_dot_array = max_dot_array
        self.method = method
        self.n_trials = n_trials
        self.save_images = save_images
        self.folder = folder
        self.area_colour = area_colour
        self.convex_hull_colour = convex_hull_colour
        self.antialiasing = antialiasing
        self.background_colour = background_colour
        self.auto_delete_image_files = auto_delete_image_files
        self.sqeeze_factor = sqeeze_factor
        self.min_numerosity = min_numerosity

    def run(self):
        cnt = 0
        da_sequence = DASequence()

        while cnt<self.n_trials:
            cnt += 1
            if da_sequence.make_by_incrementing(max_dot_array=self.max_dot_array,
                                                method=self.method,
                                                sqeeze_factor=self.sqeeze_factor,
                                                min_numerosity=self.min_numerosity):
                break
            print("remix")

        write_pil_images_of_da_sequence( dot_array_sequence=da_sequence,
                          area_colour=self.area_colour,
                          convex_hull_colour=self.convex_hull_colour,
                          antialiasing=self.antialiasing,
                          background_colour=self.background_colour)

        if self.save_images:
            pil_images_save(da_sequence_with_images=da_sequence,
                            folder= self.folder,
                            replace_images_by_filename=True)
            da_sequence.auto_delete_image_files = self.auto_delete_image_files
        else:
            da_sequence.auto_delete_image_files = False

        self.data_available.set()
        self._data_queue.put(da_sequence)