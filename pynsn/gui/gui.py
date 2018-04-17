#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

import sys
from PyQt4 import QtGui
from PIL.ImageQt import ImageQt
from .._lib.generator import DotArrayGenerator, GeneratorLogger
from .._lib.colour import Colour
from .. import pil_image


class MainWindow(QtGui.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        exitAction = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)

        #self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(exitAction)

        self.form_widget = FormWidget(self)
        self.setCentralWidget(self.form_widget)

        self.move(300, -300)
        self.setWindowTitle('xxx')
        self.form_widget.show_pict()
        self.show()


class LabeledInput(object):

    def __init__(self, label, text, width_label=180, width_edit=70):

        self.label = QtGui.QLabel(label)
        self.label.setFixedWidth(width_label)
        self.edit = QtGui.QLineEdit()
        self.edit.setFixedWidth(width_edit)
        self.text = text

    @property
    def text(self):
        return str(self.edit.text())

    @text.setter
    def text(self, v):
        self.edit.setText(str(v))

    def layout(self, vertical=False):
        if vertical:
            layout = QtGui.QVBoxLayout()
        else:
            layout = QtGui.QHBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.edit)
        return layout


class LabeledNumberInput(LabeledInput):

    def __init__(self, label, value, width_label=180, width_edit=70,
                 integer_only=True):

        LabeledInput.__init__(self, label=label, text="", width_label=width_label,
                              width_edit=width_edit)
        self.integer_only = integer_only
        if integer_only is not None:
            self.edit.setValidator(QtGui.QIntValidator())
        else:
            self.edit.setValidator(QtGui.QDoubleValidator())
        self.value = value

    @property
    def value(self):
        if self.integer_only:
            return int(self.edit.text())
        else:
            return float(self.edit.text())

    @value.setter
    def value(self, v):
        if self.integer_only:
            v = int(v)
        else:
            v = float(v)
        self.edit.setText(str(v))


class FormWidget(QtGui.QWidget):

    def __init__(self, parent):

        super(FormWidget, self).__init__(parent)

        self.logger = GeneratorLogger(log_filename="log/gui",
                                      override_log_files=True,
                                      log_colours=False,
                                      properties_different_colour=False)
        self.initUI()

    def make_stimulus(self, number):
        try:
            dot_colour = Colour(self.dot_colour.text)
        except:
            dot_colour = "green"
            self.dot_colour.text = "green"

        generator = DotArrayGenerator(
            max_array_radius=self.max_array_radius.value,
            dot_diameter_mean=self.dot_diameter_mean.value,
            dot_diameter_range=(self.dot_diameter_range_from.value,
                                self.dot_diameter_range_to.value),
            dot_diameter_std=self.dot_diameter_std.value,
            dot_colour=dot_colour,
            minimum_gap=self.minimum_gap.value,
            logger=self.logger)
        return generator.make(n_dots=number)

    def initUI(self):
        button1 = QtGui.QPushButton(" Display")
        button1.clicked.connect(self.show_pict)

        self.number = LabeledNumberInput("Number", 80)
        self.max_array_radius = LabeledNumberInput("max array radius", 200)
        self.dot_diameter_mean = LabeledNumberInput("mean diameter", 25)
        self.dot_diameter_range_from = LabeledNumberInput("diameter range from", 5)
        self.dot_diameter_range_to = LabeledNumberInput("diameter range to", 40)
        self.dot_diameter_std = LabeledNumberInput("diameter range std", 2)
        self.minimum_gap = LabeledNumberInput("minimum gap", 2)

        self.dot_colour = LabeledInput("colour", text="skyblue")
        self.colour_area = LabeledInput("Area colour", text="None")
        self.colour_convex_hull_positions = LabeledInput("Convex hull colour positions", text="None")
        self.colour_convex_hull_dots = LabeledInput("Convex hull colour dots", text="None")
        self.colour_background = LabeledInput("Background colour", text="black")
        antialiasing = True #Todo

        crtl = QtGui.QVBoxLayout()
        crtl.addWidget(button1)
        crtl.addLayout(self.number.layout())
        crtl.addLayout(self.max_array_radius.layout())
        crtl.addLayout(self.dot_diameter_mean.layout())
        crtl.addLayout(self.dot_diameter_range_from.layout())
        crtl.addLayout(self.dot_diameter_range_to.layout())
        crtl.addLayout(self.dot_diameter_std.layout())
        crtl.addLayout(self.dot_colour.layout())
        crtl.addLayout(self.minimum_gap.layout())

        crtl.addLayout(self.colour_area.layout())
        crtl.addLayout(self.colour_convex_hull_positions.layout())
        crtl.addLayout(self.colour_convex_hull_dots.layout())
        crtl.addLayout(self.colour_background.layout())


        self.picture_field = QtGui.QLabel(self)
        self.picture_field.setFixedSize(600, 600)

        hlayout = QtGui.QHBoxLayout()
        hlayout.addLayout(crtl)
        hlayout.addWidget(self.picture_field)

        self.setLayout(hlayout)

        #self.move(300, 200)
        self.setWindowTitle('Red Rock')
        self.show()

    def show_pict(self):
        da = self.make_stimulus(number=self.number.value)
        self.picture_field.setFixedSize(self.max_array_radius.value*2,
                                        self.max_array_radius.value*2)

        try:
            colour_area = Colour(self.colour_area.text)
        except:
            colour_area = None
            self.colour_area.text = "None"
        try:
            colour_convex_hull_positions = Colour(self.colour_convex_hull_positions.text)
        except:
            colour_convex_hull_positions = None
            self.colour_convex_hull_positions.text = "None"
        try:
            colour_convex_hull_dots = Colour(self.colour_convex_hull_dots.text)
        except:
            colour_convex_hull_dots = None
            self.colour_convex_hull_dots.text = "None"
        try:
            colour_background = Colour(self.colour_background.text)
        except:
            colour_background = None
            self.colour_background.text = "None"

        im = pil_image.create(da,
                         colour_area=colour_area,
                         colour_convex_hull_positions=colour_convex_hull_positions,
                         colour_convex_hull_dots=colour_convex_hull_dots,
                         colour_center_of_mass = None,
                         colour_center_of_outer_positions=None,
                         antialiasing=True,
                         colour_background=colour_background)

        qtim = ImageQt(im)
        pixmap = QtGui.QPixmap.fromImage(qtim)
        self.picture_field.setPixmap(pixmap)
        #print(dir(self.pixmap))
        self.adjustSize()
        self.parent().adjustSize()


def start():
    app = QtGui.QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    start()