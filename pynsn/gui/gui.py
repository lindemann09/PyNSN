#!/usr/bin/python

"""
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

import os
import sys
from PyQt4 import QtGui
from PIL.ImageQt import ImageQt
from .._lib.generator import DotArrayGenerator, GeneratorLogger, DASequenceGenerator
from .._lib.colour import Colour
from .. import pil_image
from .main_widget import MainWidget
from . import dialogs
from .sequence_display import SequenceDisplay

DEFAULT_ARRAY = pil_image.RandomDotImageParameter(number=40,
                                                  max_array_radius=200,
                                                  dot_colour="lime",
                                                  dot_diameter_mean=15,
                                                  dot_diameter_range=[5, 40],
                                                  dot_diameter_std=8,
                                                  minimum_gap=2,
                                                  colour_area="#3e3e3e",
                                                  colour_convex_hull_positions=None,
                                                  colour_convex_hull_dots=None,
                                                  colour_center_of_mass=None,
                                                  colour_center_of_outer_positions=None,
                                                  antialiasing=True,
                                                  colour_background="gray")

ICON = pil_image.RandomDotImageParameter(number=11,
                                         max_array_radius=200,
                                         dot_colour="lime",
                                         dot_diameter_mean=35,
                                         dot_diameter_range=[5, 80],
                                         dot_diameter_std=20,
                                         minimum_gap=2,
                                         colour_area="#3e3e3e",
                                         colour_convex_hull_positions=None,
                                         colour_convex_hull_dots="expyriment_orange",
                                         colour_center_of_mass=None,
                                         colour_center_of_outer_positions=None,
                                         antialiasing=True,
                                         colour_background=None)


class PyNSN_GUI(QtGui.QMainWindow):

    def __init__(self):

        super(PyNSN_GUI, self).__init__()

        self._image = None
        self.data_array = None
        self.set_loging(False)  # todo checkbox in settings
        self.settings = dialogs.SettingsDialog(self, DEFAULT_ARRAY)
        self.initUI()
        self.show()

    def set_loging(self, onoff):

        if onoff:
            self.logger = GeneratorLogger(log_filename="log/gui",
                                          override_log_files=True,
                                          log_colours=False,
                                          properties_different_colour=False)
        else:
            self.logger = None

    def initUI(self):

        # menues
        exitAction = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)

        saveAction = QtGui.QAction('&Save stimulus', self)
        saveAction.setShortcut('Ctrl+S')
        saveAction.triggered.connect(self.save_pixmap)

        settingsAction = QtGui.QAction('&Settings', self)
        settingsAction.triggered.connect(self.action_settings)

        printxyAction = QtGui.QAction('&Print array', self)
        printxyAction.triggered.connect(self.action_print_xy)

        matchAction = QtGui.QAction('&Match property', self)
        matchAction.triggered.connect(self.action_match)

        sequenceAction = QtGui.QAction('&Make sequence', self)
        sequenceAction.triggered.connect(self.action_make_sequence)

        # self.statusBar()
        menubar = self.menuBar()

        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(saveAction)
        fileMenu.addSeparator()
        fileMenu.addAction(exitAction)

        toolMenu = menubar.addMenu('&Tools')
        toolMenu.addAction(sequenceAction)
        toolMenu.addAction(settingsAction)
        toolMenu.addAction(matchAction)
        toolMenu.addAction(printxyAction)

        # main widget
        self.main_widget = MainWidget(self, self.settings, DEFAULT_ARRAY)
        self.setCentralWidget(self.main_widget)
        self.main_widget.btn_generate.clicked.connect(self.action_generate_btn)
        self.main_widget.dot_colour.edit.editingFinished.connect(self.action_dot_colour_change)
        self.main_widget.slider.sliderReleased.connect(self.action_slider_released)

        self.move(300, -300)
        self.setWindowTitle('PyNSN GUI')

        # ICON
        self._image, _ = pil_image.generate_random_dot_array_image(ICON, logger=None)
        self.setWindowIcon(QtGui.QIcon(self.pixmap()))
        self._image = None

        self.action_generate_btn()

    def make_new_array(self):

        para = self.get_parameter()
        generator = DotArrayGenerator(
            max_array_radius=para.max_array_radius,
            dot_diameter_mean=para.dot_diameter_mean,
            dot_diameter_range=para.dot_diameter_range,
            dot_diameter_std=para.dot_diameter_std,
            dot_colour=para.dot_colour,
            minimum_gap=para.minimum_gap,
            logger=self.logger)

        if self.settings.bicoloured.isChecked():
            n2 = self.main_widget.number2.value
            self.data_array = generator.make(n_dots=para.number + n2)
            self.data_array.features.change(indices=list(range(n2)),
                                            colour=self.main_widget.dot_colour2.text)
        else:
            self.data_array = generator.make(n_dots=para.number)

        self._image = None

    def image(self):
        if self._image is not None:
            return self._image
        else:
            para = self.get_parameter()
            self._image = pil_image.create(self.data_array,
                                           colour_area=para.colour_area,
                                           colour_convex_hull_positions=para.colour_convex_hull_positions,
                                           colour_convex_hull_dots=para.colour_convex_hull_dots,
                                           colour_center_of_mass=para.colour_center_of_mass,
                                           colour_center_of_outer_positions=para.colour_center_of_outer_positions,
                                           antialiasing=para.antialiasing,
                                           colour_background=para.colour_background)
            return self._image

    def get_parameter(self):
        # check colour input

        try:
            colour_dot = Colour(self.main_widget.dot_colour.text)
        except:
            colour_dot = DEFAULT_ARRAY.dot_colour
            self.main_widget.dot_colour.text = colour_dot
        try:
            colour_area = Colour(self.settings.colour_area.text)
        except:
            colour_area = None
            self.settings.colour_area.text = "None"
        try:
            colour_convex_hull_positions = Colour(self.settings.colour_convex_hull_positions.text)
        except:
            colour_convex_hull_positions = None
            self.settings.colour_convex_hull_positions.text = "None"
        try:
            colour_convex_hull_dots = Colour(self.settings.colour_convex_hull_dots.text)
        except:
            colour_convex_hull_dots = None
            self.settings.colour_convex_hull_dots.text = "None"
        try:
            colour_background = Colour(self.settings.colour_background.text)
        except:
            colour_background = None
            self.settings.colour_background.text = "None"

        return pil_image.RandomDotImageParameter(
            number=self.main_widget.number.value,
            max_array_radius=self.main_widget.max_array_radius.value,
            dot_colour=colour_dot,
            dot_diameter_mean=self.main_widget.dot_diameter_mean.value,
            dot_diameter_range=[self.main_widget.dot_diameter_range.value1,
                                self.main_widget.dot_diameter_range.value2],
            dot_diameter_std=self.main_widget.dot_diameter_std.value,
            minimum_gap=self.main_widget.minimum_gap.value,
            colour_area=colour_area,
            colour_convex_hull_positions=colour_convex_hull_positions,
            colour_convex_hull_dots=colour_convex_hull_dots,
            colour_center_of_mass=None,
            colour_center_of_outer_positions=None,
            antialiasing=self.settings.antialiasing.isChecked(),
            colour_background=colour_background)

    def pixmap(self):
        return QtGui.QPixmap.fromImage(ImageQt(self.image()))

    def show_current_image(self, remake_image=False):
        """"""
        if remake_image:
            self._image = None
        w = self.get_parameter().max_array_radius * 2
        self.main_widget.resize_fields(width=w, text_height=150)
        self.main_widget.picture_field.setPixmap(self.pixmap())
        self.main_widget.adjustSize()
        self.adjustSize()

    def action_generate_btn(self):
        """"""
        self.make_new_array()
        self.show_current_image()
        self.main_widget.updateUI()
        self.write_properties()

    def write_properties(self, clear_field=True):
        prop = self.data_array.get_properties()
        txt = prop.get_nice_text()
        if self.settings.bicoloured.isChecked():
            prop = self.data_array.get_properties_split_by_colours()
            for p in prop.split:
                txt += p.get_nice_text()
        if clear_field:
            self.main_widget.text_field.clear()
        self.main_widget.text_field.append(txt)

    def action_print_xy(self):
        """"""
        txt = self.current_data_array.get_csv(object_id_column=False, num_idx_column=False, colour_column=True)
        self.main_widget.text_field.append(txt)

    def save_pixmap(self):
        """"""
        # name = QtGui.QFileDialog.getSaveFileName(self, 'Save File')
        filename, extension = QtGui.QFileDialog.getSaveFileNameAndFilter(
            self, 'Save file', filter=self.tr(".png"))  # TODO multiple file formats FIXME formats selection
        self.current_image.save(filename + extension, format=str(extension[1:]).upper())

    def action_match(self):
        """"""
        prop = self.current_data_array.get_properties()
        print(dialogs.MatchPropertyDialog.get_response(self, prop))  # TODO

    def action_settings(self):
        """"""
        result = self.settings.exec_()
        self.main_widget.updateUI()
        self.show_current_image(remake_image=True)

    def action_dot_colour_change(self):
        """"""
        self.data_array.features.change(colour=self.get_parameter().dot_colour)
        self.show_current_image(remake_image=True)

    def action_slider_released(self):
        """"""
        change = self.main_widget.number.value - self.data_array.prop_numerosity
        self.data_array = self.data_array.number_deviant(change)
        self.show_current_image(remake_image=True)
        self.write_properties()
        # todo slider does not work correctly for multi colour arrays

    def action_make_sequence(self):
        match_methods, match_range = dialogs.SequenceDialog.get_response(self)
        if match_methods is not None:
            print("processing")  # Todo processing dialog
            gen = DASequenceGenerator(reference_dot_array=self.data_array, logger=self.logger)
            sequence = gen.make(match_properties=match_methods,
                                min_max_numerosity=match_range,
                                extra_space=100  # todo dialog field
                                )
            SequenceDisplay(self, da_sequence=sequence,
                            start_numerosity=self.data_array.prop_numerosity,
                            image_parameter=self.get_parameter()).exec_()


def start():
    app = QtGui.QApplication(sys.argv)
    ex = PyNSN_GUI()
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    start()
