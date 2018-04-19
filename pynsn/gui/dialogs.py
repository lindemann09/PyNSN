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


class SettingsDialog(QtGui.QDialog):

    def __init__(self, parent, default_array):

        super(SettingsDialog, self).__init__(parent)

        self.setWindowTitle("Match Dot Array Property")

        self.colour_area = misc.LabeledInput("Area", text=default_array.colour_area, case_sensitive=False)
        self.colour_background = misc.LabeledInput("Background", text=default_array.colour_background, case_sensitive=False)
        self.colour_convex_hull_positions = misc.LabeledInput("Colour positions CH", text=default_array.colour_convex_hull_positions, case_sensitive=False)
        self.colour_convex_hull_dots = misc.LabeledInput("Colour dots CH", text=default_array.colour_convex_hull_dots, case_sensitive=False)
        self.antialiasing = QtGui.QCheckBox("Antialiasing")
        self.antialiasing.setChecked(default_array.antialiasing)

        self.bicoloured = QtGui.QCheckBox("bicoloured")
        self.bicoloured.setChecked(False)


        vlayout = QtGui.QVBoxLayout()
        vlayout.addWidget(misc.heading("Colour"))
        vlayout.addLayout(self.colour_area.layout())
        vlayout.addLayout(self.colour_background.layout())

        vlayout.addWidget(misc.heading("Convex hull"))
        vlayout.addLayout(self.colour_convex_hull_positions.layout())
        vlayout.addLayout(self.colour_convex_hull_dots.layout())
        vlayout.addSpacing(10)
        vlayout.addWidget(self.antialiasing)
        vlayout.addWidget(self.bicoloured)
        vlayout.addStretch(1)


        # OK and Cancel buttons
        buttons = QtGui.QDialogButtonBox( QtGui.QDialogButtonBox.Ok, QtCore.Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)

        vlayout.addWidget(buttons)
        self.setLayout(vlayout)