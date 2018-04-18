#!/usr/bin/python

"""
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

import os
import sys
from PyQt4 import QtGui
from PIL.ImageQt import ImageQt
from .._lib.generator import DotArrayGenerator, GeneratorLogger
from .. import pil_image
from .main_widget import MainWidget
from .dialogs import MatchPropertyDialog
from .. import pil_image


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

        self.set_loging(False) #False
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

        printxyAction= QtGui.QAction('&Print array', self)
        printxyAction.triggered.connect(self.action_print_xy)

        matchAction= QtGui.QAction('&Match property', self)
        matchAction.triggered.connect(self.action_match)

        #self.statusBar()
        menubar = self.menuBar()

        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(saveAction)
        fileMenu.addSeparator()
        fileMenu.addAction(exitAction)

        toolMenu = menubar.addMenu('&Tools')
        toolMenu.addAction(matchAction)
        toolMenu.addAction(printxyAction)

        #main widget
        self.main_widget = MainWidget(self)
        self.setCentralWidget(self.main_widget)

        #actions
        self.main_widget.btn_display.clicked.connect(self.action_display_btn)

        self.move(300, -300)
        self.setWindowTitle('PyNSN GUI')

        # dot array
        pixmap = self.make_image(ICON)
        self.setWindowIcon(QtGui.QIcon(pixmap))

        self.action_display_btn()


    def make_image(self, para):
        """returns pix map"""
        self.current_image, self.current_data_array = pil_image.generate_random_dot_array_image(para, logger=self.logger)
        return QtGui.QPixmap.fromImage(ImageQt(self.current_image))

    def show_pixmap(self, pixmap):
        self.main_widget.picture_field.setPixmap(pixmap)
        self.main_widget.adjustSize()
        self.adjustSize()

    def action_display_btn(self):
        para = self.main_widget.all_parameter
        pixmap = self.make_image(para)
        self.main_widget.resize_fields(width= para.max_array_radius*2,
                                       text_height=150)
        self.show_pixmap(pixmap)

        prop = self.current_data_array.get_properties()
        txt = prop.get_nice_text()
        self.main_widget.text_field.append(txt)

    def action_print_xy(self):

        txt = self.current_data_array.get_csv(object_id_column=False, num_idx_column=False, colour_column=True)
        self.main_widget.text_field.append(txt)


    def save_pixmap(self):
        #name = QtGui.QFileDialog.getSaveFileName(self, 'Save File')
        filename, extension = QtGui.QFileDialog.getSaveFileNameAndFilter(
            self, 'Save file', filter=self.tr(".png")) #TODO multiple file formats FIXME formats selection
        print(extension)
        self.current_image.save(filename, format=extension[1:].upper())

    def action_match(self):
        prop = self.current_data_array.get_properties()
        print(MatchPropertyDialog.get_response(self, prop))

def start():
    app = QtGui.QApplication(sys.argv)
    ex = PyNSN_GUI()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    start()