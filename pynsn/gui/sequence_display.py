"""
"""

from __future__ import absolute_import
from builtins import zip, filter, range, super

from PyQt4 import QtGui,  QtCore
from PIL.ImageQt import ImageQt
from .. import pil_image
from . import misc

class SequenceDisplay(QtGui.QDialog):

    def __init__(self, parent, da_sequence, start_numerosity, image_parameter):

        super(SequenceDisplay, self).__init__(parent)

        self.setWindowTitle("Dot Array Sequence")
        self.da_sequence = da_sequence
        self.pixmap_width = da_sequence.dot_arrays[0].max_array_radius * 2

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
        self.pixmaps = []
        for da in misc.progressbar(self.da_sequence.dot_arrays, label="make images",
                              win_title="Dot Array Sequence"): #todo cancle iterator not tested
            im = pil_image.create(da,
                                           colour_area=image_parameter.colour_area,
                                           colour_convex_hull_positions=image_parameter.colour_convex_hull_positions,
                                           colour_convex_hull_dots=image_parameter.colour_convex_hull_dots,
                                           colour_center_of_mass=image_parameter.colour_center_of_mass,
                                           colour_center_of_outer_positions=image_parameter.colour_center_of_outer_positions,
                                           antialiasing=image_parameter.antialiasing,
                                           colour_background=image_parameter.colour_background)
            self.pixmaps.append(QtGui.QPixmap.fromImage(ImageQt(im)))

        self.updateUI()

    def updateUI(self):
        num = self.slider.value()
        #print("{}-{}".format(num, self.da_sequence.dot_arrays[self.da_sequence.numerosity_idx[num]].prop_numerosity)) #todo depict somewhere
        im = self.pixmaps[self.da_sequence.numerosity_idx[num]]

        self.picture_field.setPixmap(im)
        self.adjustSize()


    def action_slider_change(self):
        self.updateUI()


