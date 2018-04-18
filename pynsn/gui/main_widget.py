"""
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

from PyQt4 import QtGui
from .. import pil_image
from .misc import heading, LabeledNumberInput, LabeledNumberInputTwoValues, LabeledInput
from .._lib.colour import Colour

DEFAULT_ARRAY = pil_image.RandomDotImageParameter(number=40,
                           max_array_radius=200,
                           dot_colour="lime",
                           dot_diameter_mean=25,
                           dot_diameter_range=[5, 40],
                           dot_diameter_std=8,
                           minimum_gap=2,
                           colour_area="#3e3e3e",
                           colour_convex_hull_positions=None,
                           colour_convex_hull_dots="skyblue",
                           colour_center_of_mass=None,
                           colour_center_of_outer_positions=None,
                           antialiasing=True,
                           colour_background="gray")

class MainWidget(QtGui.QWidget):

    def __init__(self, parent):

        super(MainWidget, self).__init__(parent)
        self.initUI()


    def initUI(self):

        self.btn_display = QtGui.QPushButton(" Display")

        self.number = LabeledNumberInput("Number", DEFAULT_ARRAY.number)
        self.max_array_radius = LabeledNumberInput("Max radius", DEFAULT_ARRAY.max_array_radius)
        self.dot_diameter_mean = LabeledNumberInput("Mean diameter", DEFAULT_ARRAY.dot_diameter_mean)
        self.dot_diameter_std = LabeledNumberInput("Diameter range std", DEFAULT_ARRAY.dot_diameter_std)
        self.dot_diameter_range = LabeledNumberInputTwoValues("Diameter range from",
                                                              value1=DEFAULT_ARRAY.dot_diameter_range[0],
                                                              value2=DEFAULT_ARRAY.dot_diameter_range[1])

        self.minimum_gap = LabeledNumberInput("Minimum gap", DEFAULT_ARRAY.minimum_gap)

        self.dot_colour = LabeledInput("Colour", text=DEFAULT_ARRAY.dot_colour, case_sensitive=False)
        self.colour_area = LabeledInput("Area", text=DEFAULT_ARRAY.colour_area, case_sensitive=False)
        self.colour_background = LabeledInput("Background", text=DEFAULT_ARRAY.colour_background, case_sensitive=False)
        self.colour_convex_hull_positions = LabeledInput("Colour positions CH", text=DEFAULT_ARRAY.colour_convex_hull_positions, case_sensitive=False)
        self.colour_convex_hull_dots = LabeledInput("Colour dots CH", text=DEFAULT_ARRAY.colour_convex_hull_dots, case_sensitive=False)
        antialiasing = DEFAULT_ARRAY.colour_convex_hull_dots #Todo

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

        # Add text field
        fields = QtGui.QVBoxLayout()
        self.text_field = QtGui.QTextBrowser(self)
        self.text_field.setFont(QtGui.QFont("Courier New", 8))
        self.text_field.setReadOnly(True)
        self.text_field.setWordWrapMode(False)
        self.text_field.setVerticalScrollBar(QtGui.QScrollBar())
        self.text_field.setHorizontalScrollBar(QtGui.QScrollBar())
        self.picture_field = QtGui.QLabel(self)
        self.resize_fields(400)
        fields.addWidget(self.picture_field)
        fields.addWidget(self.text_field)


        hlayout = QtGui.QHBoxLayout()
        hlayout.addLayout(ctrl)
        hlayout.addLayout(fields)
        self.setLayout(hlayout)

    def resize_fields(self, width, text_height=150, minium_text_width=300):
        self.picture_field.setFixedSize(width, width)
        if width<minium_text_width:
            width=minium_text_width
        self.text_field.setFixedSize(width, text_height)


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

        return pil_image.RandomDotImageParameter(
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