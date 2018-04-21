"""
"""

from __future__ import absolute_import
from builtins import zip, filter, range, super

from PyQt4 import QtGui,  QtCore
from . import misc
from .._lib import continuous_property as cp


class SequenceDisplay(QtGui.QDialog):

    def __init__(self, parent, da_sequence):

        super(SequenceDisplay, self).__init__(parent)

        self.setWindowTitle("Dot Array Sequence")
        self.da_sequence = da_sequence

        self.picture_field = QtGui.QLabel(self)

        self.slider = QtGui.QSlider(QtCore.Qt.Vertical)
        num_range = da_sequence.min_max_numerosity
        self.slider.setMinimum(num_range[0])
        self.slider.setMaximum(num_range[1])
        self.slider.setValue(da_sequence.prop_numerosity)
        self.slider.valueChanged.connect(self.action_slider_change)

        hlayout = QtGui.QHBoxLayout()
        hlayout.addLayout(self.picture_field)
        hlayout.addWidget(self.slider)
        self.setLayout(hlayout)
        self.updateUI()

    def updateUI(self):
        pass

    def action_slider_change(self):
        pass


