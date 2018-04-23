"""
"""

from __future__ import absolute_import
from builtins import zip, filter, range, super

from PyQt4 import QtGui, QtCore
from . import misc
from .._lib import continuous_property as cp


class MatchPropertyDialog(QtGui.QDialog):
    def __init__(self, parent, properties):
        super(MatchPropertyDialog, self).__init__(parent)

        self.setWindowTitle("Match Dot Array Property")

        self.properties = properties
        self.comboBox = QtGui.QComboBox(self)
        self.comboBox.addItem(cp.DotDiameter().long_label)  # 0
        self.comboBox.addItem(cp.Density().long_label)  # 1
        self.comboBox.addItem(cp.ConvexHull().long_label)  # 2
        self.comboBox.addItem(cp.SurfaceArea().long_label)  # 3
        self.comboBox.addItem(cp.Circumference().long_label)  # 4
        self.comboBox.activated[str].connect(self.choice)

        self.num_input = misc.NumberInput(width_edit=150, value=0)
        self.choice(cp.DotDiameter().long_label)

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

        if selection == cp.DotDiameter().long_label:
            self.num_input.value = self.properties.mean_dot_diameter
        elif selection == cp.Density().long_label:
            self.num_input.value = self.properties.density
        elif selection == cp.ConvexHull().long_label:
            self.num_input.value = self.properties.convex_hull_area
        elif selection == cp.Circumference().long_label:
            self.num_input.value = self.properties.total_circumference
        elif selection == cp.SurfaceArea().long_label:
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
    def __init__(self, parent, image_parameter):
        super(SettingsDialog, self).__init__(parent)

        self.setWindowTitle("Dot Array Property")

        self.colour_area = misc.LabeledInput("Area",
                                             text=image_parameter.colour_area.colour,
                                             case_sensitive=False)
        self.colour_background = misc.LabeledInput("Background",
                                                   text=image_parameter.colour_background.colour,
                                                   case_sensitive=False)
        self.colour_convex_hull_positions = misc.LabeledInput("Colour positions CH",
                                                              text=image_parameter.colour_convex_hull_positions.colour,
                                                              case_sensitive=False)
        self.colour_convex_hull_dots = misc.LabeledInput("Colour dots CH",
                                                         text=image_parameter.colour_convex_hull_dots.colour,
                                                         case_sensitive=False)
        self.antialiasing = QtGui.QCheckBox("Antialiasing")
        self.antialiasing.setChecked(image_parameter.antialiasing)

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
        buttons = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok, QtCore.Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)

        vlayout.addWidget(buttons)
        self.setLayout(vlayout)


class SequenceDialog(QtGui.QDialog):
    def __init__(self, parent):

        super(SequenceDialog, self).__init__(parent)

        self.match_methods = []

        self.setWindowTitle("Sequence Dialog")

        self.match_diameter = QtGui.QCheckBox(cp.DotDiameter.long_label)
        self.match_area = QtGui.QCheckBox(cp.SurfaceArea.long_label)
        self.match_circumference = QtGui.QCheckBox(cp.Circumference.long_label)
        self.match_density = QtGui.QCheckBox(cp.Density.long_label)
        self.match_convex_hull = QtGui.QCheckBox(cp.ConvexHull.long_label)
        self.match_ch_presision = misc.LabeledNumberInput("Convex_hull presision",
                                                          value=cp.ConvexHull().match_presision,
                                                          integer_only=False)
        self.match_density_ratio = misc.LabeledNumberInput("Ratio convex_hull/area",
                                                           value=cp.Density().match_ratio_convhull2area,
                                                           integer_only=False, min=0, max=1)
        self.match_range = misc.LabeledNumberInputTwoValues("Sequence Range", value1=10, value2=100)
        self.match_extra_space = misc.LabeledNumberInput("Extra space",
                                                           value=50, integer_only=True, min=0)

        self.match_area.toggled.connect(self.ui_update)
        self.match_convex_hull.toggled.connect(self.ui_update)
        self.match_diameter.toggled.connect(self.ui_update)
        self.match_circumference.toggled.connect(self.ui_update)
        self.match_density.toggled.connect(self.ui_update)
        self.match_ch_presision.edit.editingFinished.connect(self.ui_update)
        self.match_density_ratio.edit.editingFinished.connect(self.ui_update)

        # OK and Cancel buttons
        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vlayout = QtGui.QVBoxLayout()
        vlayout.addWidget(misc.heading("Matching"))
        vlayout.addWidget(self.match_diameter)
        vlayout.addWidget(self.match_area)
        vlayout.addWidget(self.match_circumference)
        vlayout.addWidget(self.match_convex_hull)
        vlayout.addWidget(self.match_density)
        vlayout.addSpacing(10)
        vlayout.addWidget(misc.heading("Matching parameter"))
        vlayout.addLayout(self.match_ch_presision.layout())
        vlayout.addLayout(self.match_density_ratio.layout())
        vlayout.addLayout(self.match_range.layout())
        vlayout.addSpacing(10)
        vlayout.addLayout(self.match_extra_space.layout())
        vlayout.addSpacing(20)

        vlayout.addWidget(buttons)
        self.setLayout(vlayout)
        self.ui_update()

    def ui_update(self):
        # get methods
        selected = []
        all = [cp.DotDiameter()]
        if self.match_diameter.isChecked():
            selected.append(all[-1])

        all.append(cp.Circumference())
        if self.match_circumference.isChecked():
            selected.append(all[-1])

        all.append(cp.SurfaceArea())
        if self.match_area.isChecked():
            selected.append(all[-1])

        all.append(cp.ConvexHull(match_presision=self.match_ch_presision.value))
        if self.match_convex_hull.isChecked():
            selected.append(all[-1])

        all.append(cp.Density(match_ratio_convhull2area=self.match_density_ratio.value,
                              convex_hull_precision=self.match_ch_presision.value))
        if self.match_density.isChecked():
            selected.append(all[-1])

        self.match_diameter.setEnabled(True)
        self.match_area.setEnabled(True)
        self.match_circumference.setEnabled(True)
        self.match_density.setEnabled(True)
        self.match_convex_hull.setEnabled(True)
        for x in all:
            if x not in selected:
                if sum(map(lambda s: s.is_dependent(x), selected)) > 0:  # any dependency
                    if isinstance(x, cp.DotDiameter):
                        self.match_diameter.setEnabled(False)
                        self.match_diameter.setChecked(False)
                    if isinstance(x, cp.SurfaceArea):
                        self.match_area.setEnabled(False)
                        self.match_area.setChecked(False)
                    if isinstance(x, cp.Circumference):
                        self.match_circumference.setEnabled(False)
                        self.match_circumference.setChecked(False)
                    if isinstance(x, cp.ConvexHull):
                        self.match_convex_hull.setEnabled(False)
                        self.match_convex_hull.setChecked(False)
                    if isinstance(x, cp.Density):
                        self.match_density.setEnabled(False)
                        self.match_density.setChecked(False)

        self.match_methods = selected

    @staticmethod
    def get_response(parent):

        dialog = SequenceDialog(parent)
        result = dialog.exec_()
        if result == QtGui.QDialog.Accepted:

            return (dialog.match_methods, [dialog.match_range.value1, dialog.match_range.value2],
                    dialog.match_extra_space.value)
        else:
            return (None, None, None)
