"""
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

from PyQt4 import QtGui,  QtCore
from .parameter import RandomDotArrayImageAllParameter
from .._lib.colour import Colour

DEFAULT = RandomDotArrayImageAllParameter(number=40,
                           max_array_radius=200,
                           dot_colour="skyblue",
                           dot_diameter_mean=25,
                           dot_diameter_range=[5, 40],
                           dot_diameter_std=8,
                           minimum_gap=2,
                           colour_area="#3e3e3e",
                           colour_convex_hull_positions=None,
                           colour_convex_hull_dots="yellow",
                           colour_center_of_mass=None,
                           colour_center_of_outer_positions=None,
                           antialiasing=True,
                           colour_background="gray")

class LabeledInput(object):

    def __init__(self, label, text, width_label=180, width_edit=70, case_sensitive=True):

        self.label = QtGui.QLabel(label)
        self.label.setFixedWidth(width_label)
        self.edit = QtGui.QLineEdit()
        self.edit.setFixedWidth(width_edit)
        self.edit.setAlignment(QtCore.Qt.AlignRight)
        self.case_sensitive = case_sensitive
        self.text = text

    @property
    def text(self):
        rtn = str(self.edit.text())
        if not self.case_sensitive:
            return rtn.lower()
        else:
            return rtn

    @text.setter
    def text(self, v):
        if not self.case_sensitive:
            self.edit.setText(str(v).lower())
        else:
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

class LabeledNumberInputTwoValues(object):

    def __init__(self, label, value1, value2, width_label=180, width_edits=35,
                 integer_only=True):

        self.input1 = LabeledNumberInput(label=label, value=value1, width_label=width_label, width_edit=width_edits,
                                         integer_only=integer_only)

        self.input2 = LabeledNumberInput(label="", value=value2, width_label=width_label, width_edit=width_edits,
                                         integer_only=integer_only)

    @property
    def value1(self):
        return self.input1.value

    @value1.setter
    def value1(self, v):
        self.input1.value = v

    @property
    def value2(self):
        return self.input2.value

    @value2.setter
    def value2(self, v):
        self.input2.value = v

    def layout(self, vertical=False):
        layout = self.input1.layout(vertical)
        layout.addWidget(self.input2.edit)
        return layout


def heading(text):

    boldFont = QtGui.QFont()
    boldFont.setBold(True)
    heading = QtGui.QLabel(text)
    heading.setFont(boldFont)
    heading.setStyleSheet("QLabel { color : black; }")
    return heading


class MainWidget(QtGui.QWidget):

    def __init__(self, parent):

        super(MainWidget, self).__init__(parent)
        self.initUI()


    def initUI(self):

        self.btn_display = QtGui.QPushButton(" Display")

        self.number = LabeledNumberInput("Number", DEFAULT.number)
        self.max_array_radius = LabeledNumberInput("Max radius", DEFAULT.max_array_radius)
        self.dot_diameter_mean = LabeledNumberInput("Mean diameter", DEFAULT.dot_diameter_mean)
        self.dot_diameter_std = LabeledNumberInput("Diameter range std", DEFAULT.dot_diameter_std)
        self.dot_diameter_range = LabeledNumberInputTwoValues("Diameter range from",
                                                              value1=DEFAULT.dot_diameter_range[0],
                                                              value2=DEFAULT.dot_diameter_range[1])

        self.minimum_gap = LabeledNumberInput("Minimum gap", DEFAULT.minimum_gap)

        self.dot_colour = LabeledInput("Colour", text=DEFAULT.dot_colour, case_sensitive=False)
        self.colour_area = LabeledInput("Area", text=DEFAULT.colour_area, case_sensitive=False)
        self.colour_background = LabeledInput("Background", text=DEFAULT.colour_background, case_sensitive=False)
        self.colour_convex_hull_positions = LabeledInput("Colour positions CH", text=DEFAULT.colour_convex_hull_positions, case_sensitive=False)
        self.colour_convex_hull_dots = LabeledInput("Colour dots CH", text=DEFAULT.colour_convex_hull_dots, case_sensitive=False)
        antialiasing = DEFAULT.colour_convex_hull_dots #Todo

        ctrl = QtGui.QVBoxLayout()
        ctrl.addWidget(self.btn_display)
        ctrl.addSpacing(10)
        ctrl.addLayout(self.number.layout())

        ctrl.addWidget(heading("Dot"))
        ctrl.addLayout(self.dot_diameter_mean.layout())
        ctrl.addLayout(self.dot_diameter_range.layout())
        ctrl.addLayout(self.dot_diameter_std.layout())
        ctrl.addLayout(self.dot_colour.layout())

        ctrl.addWidget(heading("Array"))
        ctrl.addLayout(self.max_array_radius.layout())
        ctrl.addLayout(self.minimum_gap.layout())

        ctrl.addWidget(heading("Colour"))
        ctrl.addLayout(self.colour_area.layout())
        ctrl.addLayout(self.colour_background.layout())

        ctrl.addWidget(heading("Convex hull"))
        ctrl.addLayout(self.colour_convex_hull_positions.layout())
        ctrl.addLayout(self.colour_convex_hull_dots.layout())
        ctrl.addStretch(1)

        self.picture_field = QtGui.QLabel(self)
        self.picture_field.setFixedSize(600, 600)

        hlayout = QtGui.QHBoxLayout()
        hlayout.addLayout(ctrl)
        hlayout.addWidget(self.picture_field)

        self.setLayout(hlayout)

        self.show() #TODO required?

    @property
    def all_parameter(self):
        # check colour input
        try:
            colour_dot = Colour(self.dot_colour.text)
        except:
            colour_dot = "green"
            self.dot_colour.text = colour_dot
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

        return RandomDotArrayImageAllParameter(
                           number=self.number.value,
                           max_array_radius=self.max_array_radius.value,
                           dot_colour=colour_dot,
                           dot_diameter_mean=self.dot_diameter_mean.value,
                           dot_diameter_range=[self.dot_diameter_range.value1, self.dot_diameter_range.value2],
                           dot_diameter_std= self.dot_diameter_std.value,
                           minimum_gap= self.minimum_gap.value,
                           colour_area= colour_area,
                           colour_convex_hull_positions=colour_convex_hull_positions,
                           colour_convex_hull_dots=colour_convex_hull_dots,
                           colour_center_of_mass= None,
                           colour_center_of_outer_positions=None,
                           antialiasing=True,
                           colour_background= colour_background)