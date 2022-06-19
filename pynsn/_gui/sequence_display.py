"""
"""
from multiprocessing import Pool

from PIL.ImageQt import ImageQt
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QDialog, QHBoxLayout, QLabel, QSlider

from ..image import pil
from . import misc


def _map_make_image(x):
    da, colours, aa = x
    return pil.create(object_array=da, colours=colours,
                      antialiasing=aa)


class SequenceDisplay(QDialog):

    def __init__(self, parent, da_sequence, start_numerosity, image_colours,
                 antialiasing):
        super(SequenceDisplay, self).__init__(parent)

        self.setWindowTitle("Dot Array Sequence")
        self.da_sequence = da_sequence
        self.pixmap_width = da_sequence.dot_arrays[0].target_array_radius * 2

        self.picture_field = QLabel(self)
        self.picture_field.setFixedSize(self.pixmap_width, self.pixmap_width)

        self.slider = QSlider(Qt.Vertical)
        num_range = da_sequence.min_max_numerosity
        self.slider.setMinimum(num_range[0])
        self.slider.setMaximum(num_range[1])
        self.slider.setValue(start_numerosity)
        self.slider.valueChanged.connect(self.action_slider_change)

        hlayout = QHBoxLayout()
        hlayout.addWidget(self.picture_field)
        hlayout.addWidget(self.slider)
        self.setLayout(hlayout)

        # make images
        image_colours = [image_colours] * len(self.da_sequence.dot_arrays)
        antialiasing = [antialiasing] * len(self.da_sequence.dot_arrays)

        iter_images = Pool().imap(_map_make_image,
                                  zip(self.da_sequence.dot_arrays,
                                      image_colours, antialiasing))
        progbar_iter = misc.progressbar_iterator(iter_images,
                                                 n_elements=len(self.da_sequence.dot_arrays),
                                                 label="make images", win_title="Dot Array Sequence")

        self.pixmaps = list(map(lambda im: QPixmap.fromImage(ImageQt(im)), progbar_iter))
        self.updateUI()

    def updateUI(self):
        num = self.slider.value()
        idx = self.da_sequence.numerosity_idx[num]
        feat = self.da_sequence.dot_arrays[idx]._features.get_features_text(
            extended_format=False, with_hash=False)
        self.setWindowTitle(feat)
        self.picture_field.setPixmap(self.pixmaps[idx])
        self.adjustSize()

    def action_slider_change(self):
        self.updateUI()
