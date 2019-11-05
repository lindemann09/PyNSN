"""
"""
from multiprocessing import Pool
from PyQt4 import QtGui, QtCore
from PIL.ImageQt import ImageQt
from pynsn.lib import pil_image
from . import misc

def _map_make_image(x):
    da, colours, aa = x
    return pil_image.create(dot_array=da, colours=colours,
                            antialiasing=aa)


class SequenceDisplay(QtGui.QDialog):

    def __init__(self, parent, da_sequence, start_numerosity, image_parameter):
        super(SequenceDisplay, self).__init__(parent)

        self.setWindowTitle("Dot Array Sequence")
        self.da_sequence = da_sequence
        self.pixmap_width = da_sequence.dot_arrays[0].target_array_radius * 2

        self.picture_field = QtGui.QLabel(self)
        self.picture_field.setFixedSize(self.pixmap_width, self.pixmap_width)

        self.slider = QtGui.QSlider(QtCore.Qt.Vertical)
        num_range = da_sequence.min_max_numerosity
        self.slider.setMinimum(num_range[0])
        self.slider.setMaximum(num_range[1])
        self.slider.setValue(start_numerosity)
        self.slider.valueChanged.connect(self.action_slider_change)

        hlayout = QtGui.QHBoxLayout()
        hlayout.addWidget(self.picture_field)
        hlayout.addWidget(self.slider)
        self.setLayout(hlayout)

        # make images
        image_colours = pil_image.ImageColours(
            target_area=image_parameter.colour_target_area,
            field_area=image_parameter.colour_field_area,
            field_area_outer=image_parameter.colour_field_area_outer,
            center_of_mass=image_parameter.colour_center_of_mass,
            center_of_outer_positions=image_parameter.colour_center_of_outer_positions,
            background=image_parameter.colour_background)
        image_colours = [image_colours] * len(self.da_sequence.dot_arrays)
        antialiasing = [image_parameter.antialiasing] * len(self.da_sequence.dot_arrays)

        iter_images = Pool().imap(_map_make_image,
                                  zip(self.da_sequence.dot_arrays,
                                      image_colours, antialiasing))
        progbar_iter = misc.progressbar_iterator(iter_images,
                                                 n_elements=len(self.da_sequence.dot_arrays),
                                                 label="make images", win_title="Dot Array Sequence")

        self.pixmaps = list(map(lambda im: QtGui.QPixmap.fromImage(ImageQt(im)), progbar_iter))
        self.updateUI()

    def updateUI(self):
        num = self.slider.value()
        idx = self.da_sequence.numerosity_idx[num]
        feat = self.da_sequence.dot_arrays[idx].get_features_text(extended_format=False, with_object_id=False)
        self.setWindowTitle(feat)
        self.picture_field.setPixmap(self.pixmaps[idx])
        self.adjustSize()

    def action_slider_change(self):
        self.updateUI()
