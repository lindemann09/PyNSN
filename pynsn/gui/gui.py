#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
ZetCode PyQt5 tutorial

In this example, we position two push
buttons in the bottom-right corner
of the window.

Author: Jan Bodnar
Website: zetcode.com
Last edited: August 2017
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

import sys
import os
from PyQt4 import QtGui
from PIL.ImageQt import ImageQt
from .._lib.generator import DotArrayGenerator, GeneratorLogger
from .. import pil_image


class MainWindow(QtGui.QMainWindow):

    def __init__(self, pict_window):
        super(MainWindow, self).__init__()
        exitAction = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)

        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(exitAction)

        self.pict_window = pict_window
        self.form_widget = FormWidget(self)
        self.setCentralWidget(self.form_widget)

        self.setGeometry(300, 300, 500, 600)
        self.setWindowTitle('Menubar')
        self.show()

class FormWidgetFull(QtGui.QWidget):

    def __init__(self, parent):

        super(FormWidgetFull, self).__init__(parent)

        self.logger = GeneratorLogger(log_filename="log/gui",
                                      override_log_files=True,
                                      log_colours=False,
                                      properties_different_colour=False)
        self.initUI()

    def make_stimulus(self, number):
        generator = DotArrayGenerator(
            max_array_radius=300,
            dot_diameter_mean=25,
            dot_diameter_range=(5, 40),
            dot_diameter_std=2,
            dot_colour="skyblue",
            minimum_gap=5,
            logger=self.logger)
        return generator.make(n_dots=number)

    def initUI(self):
        button1 = QtGui.QPushButton("Button 1")
        button1.clicked.connect(self.show_pict)


        layout = QtGui.QVBoxLayout(self)
        self.pixmap = QtGui.QPixmap()

        self.pic = QtGui.QLabel(self)
        self.pic.setFixedSize(600, 600)

        layout.addWidget(button1)
        layout.addWidget(self.pic)

        self.setLayout(layout)

        #self.move(300, 200)
        self.setWindowTitle('Red Rock')
        self.show()

    def show_pict(self):
        da = self.make_stimulus(number=80)
        im = pil_image.create(da,
                         colour_area=None,
                         colour_convex_hull_positions=None,
                         colour_convex_hull_dots=None,
                         colour_center_of_mass = None,
                         colour_center_of_outer_positions=None,
                         antialiasing=True,
                         colour_background=(0, 0, 0),
                         default_dot_colour="lightgreen")


        qtim = ImageQt(im)
        pixmap = QtGui.QPixmap.fromImage(qtim)
        self.pic.setPixmap(pixmap)
        #print(dir(self.pixmap))


class FormWidget(QtGui.QWidget):

    def __init__(self, parent):

        super(FormWidget, self).__init__(parent)

        self.logger = GeneratorLogger(log_filename="log/gui",
                                      override_log_files=True,
                                      log_colours=False,
                                      properties_different_colour=False)
        self.initUI()

    def make_stimulus(self, number):
        generator = DotArrayGenerator(
            max_array_radius=300,
            dot_diameter_mean=25,
            dot_diameter_range=(5, 40),
            dot_diameter_std=2,
            dot_colour="skyblue",
            minimum_gap=5,
            logger=self.logger)
        return generator.make(n_dots=number)

    def initUI(self):
        button1 = QtGui.QPushButton("Button 1")
        button1.clicked.connect(self.show_pict)


        layout = QtGui.QVBoxLayout(self)
        layout.addWidget(button1)
        self.setLayout(layout)

        #self.move(300, 200)
        self.setWindowTitle('Red Rock')
        self.show()


    def show_pict(self):
        da = self.make_stimulus(number=80)
        im = pil_image.create(da,
                         colour_area=None,
                         colour_convex_hull_positions=None,
                         colour_convex_hull_dots=None,
                         colour_center_of_mass = None,
                         colour_center_of_outer_positions=None,
                         antialiasing=True,
                         colour_background=(0, 0, 0),
                         default_dot_colour="lightgreen")


        qtim = ImageQt(im)
        pixmap = QtGui.QPixmap.fromImage(qtim)
        #self.pic.setPixmap(pixmap)
        #print(dir(self.pixmap))



class PictureWindow(QtGui.QWidget):

    def __init__(self, parent=None):

        super(PictureWindow, self).__init__(parent)
        self.initUI()

    def initUI(self):

        self.pic = QtGui.QLabel(self)
        self.pic.setFixedSize(600, 600)

        layout = QtGui.QVBoxLayout(self)
        layout.addWidget(self.pic)
        self.setLayout(layout)

        #self.move(300, 200)
        self.setWindowTitle('XX')
        self.show()


def start():
    app = QtGui.QApplication(sys.argv)
    pw = PictureWindow()
    ex = MainWindow(pict_window=pw)
    pw.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    start()