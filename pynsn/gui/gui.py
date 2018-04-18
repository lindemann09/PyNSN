#!/usr/bin/python

"""
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

import sys
from PyQt4 import QtGui
from PIL.ImageQt import ImageQt
from .._lib.generator import DotArrayGenerator, GeneratorLogger
from .. import pil_image
from .layout import MainWidget
from .parameter import ICON


class PyNSN_GUI(QtGui.QMainWindow):

    def __init__(self):
        super(PyNSN_GUI, self).__init__()

        self.logger = GeneratorLogger(log_filename="log/gui",
                                      override_log_files=True,
                                      log_colours=False,
                                      properties_different_colour=False)
        self.initUI()
        self.show()

    def initUI(self):
        # menues
        exitAction = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)
        #self.statusBar()
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(exitAction)

        #main widget
        self.main_widget = MainWidget(self)
        self.setCentralWidget(self.main_widget)

        #actions
        self.main_widget.btn_display.clicked.connect(self.action_display_btn)

        self.move(300, -300)
        self.setWindowTitle('PyNSN GUI')

        # dot array
        pixmap = self.make_pixmap(ICON)
        self.setWindowIcon(QtGui.QIcon(pixmap))

        self.action_display_btn()


    def make_pixmap(self, para):

        generator = DotArrayGenerator(
            max_array_radius= para.max_array_radius,
            dot_diameter_mean=para.dot_diameter_mean,
            dot_diameter_range=para.dot_diameter_range,
            dot_diameter_std=para.dot_diameter_std,
            dot_colour=para.dot_colour,
            minimum_gap=para.minimum_gap,
            logger=self.logger)

        da = generator.make(n_dots=para.number)
        self.main_widget.picture_field.setFixedSize(para.max_array_radius*2,
                                        para.max_array_radius*2)
        im = pil_image.create(da,
                         colour_area=para.colour_area,
                         colour_convex_hull_positions=para.colour_convex_hull_positions,
                         colour_convex_hull_dots=para.colour_convex_hull_dots,
                         colour_center_of_mass = para.colour_center_of_mass,
                         colour_center_of_outer_positions=para.colour_center_of_outer_positions,
                         antialiasing=para.colour_center_of_outer_positions,
                         colour_background=para.colour_background)

        return QtGui.QPixmap.fromImage(ImageQt(im))

    def show_pixmap(self, pixmap):
        self.main_widget.picture_field.setPixmap(pixmap)
        self.main_widget.adjustSize()
        self.adjustSize()

    def action_display_btn(self):
        pixmap = self.make_pixmap(self.main_widget.all_parameter)
        self.show_pixmap(pixmap)



def start():
    app = QtGui.QApplication(sys.argv)
    ex = PyNSN_GUI()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    start()