"""
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

from PyQt4 import QtGui,  QtCore
from . import misc
from .._lib import constants as const


class MatchPropertyDialog(QtGui.QDialog):
    def __init__(self, parent, properties):
        super(MatchPropertyDialog, self).__init__(parent)

        self.setWindowTitle("Match Dot Array Property")

        self.properties = properties
        self.comboBox = QtGui.QComboBox(self)
        self.comboBox.addItem(const.P_DIAMETER) #0
        self.comboBox.addItem(const.P_DENSITY) #1
        self.comboBox.addItem(const.P_CONVEX_HULL) #2
        self.comboBox.addItem(const.P_AREA) #3
        self.comboBox.addItem(const.P_CIRCUMFERENCE) #4
        self.comboBox.activated[str].connect(self.choice)

        self.num_input = misc.NumberInput(width_edit=150, value=0)
        self.choice(const.P_DIAMETER)

        vlayout = QtGui.QVBoxLayout(self)
        hlayout = QtGui.QHBoxLayout()

        hlayout.addWidget(self.comboBox)
        hlayout.addWidget(self.num_input.edit)

        # OK and Cancel buttons
        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vlayout.addLayout(hlayout)
        vlayout.addWidget(buttons)

    def choice(self, selection):

        if selection == const.P_DIAMETER:
            self.num_input.value = self.properties.mean_dot_diameter
        elif selection == const.P_DENSITY:
            self.num_input.value = self.properties.density
        elif selection == const.P_CONVEX_HULL:
            self.num_input.value = self.properties.convex_hull_area
        elif selection == const.P_CIRCUMFERENCE:
            self.num_input.value = self.properties.total_circumference
        elif selection == const.P_AREA:
            self.num_input.value = self.properties.total_surface_area

    @staticmethod
    def get_response(parent, prop):

        dialog = MatchPropertyDialog(parent, prop)
        result = dialog.exec_()
        if result == QtGui.QDialog.Accepted:
            return (True, True)
        else:
            return (None, None)