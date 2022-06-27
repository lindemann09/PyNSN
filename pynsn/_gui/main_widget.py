"""
"""
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QWidget, QLabel, \
    QTextBrowser, QScrollBar, QPushButton

from .misc import heading, LabeledNumberInput, LabeledNumberInputTwoValues, \
    LabeledInput
from ..factory import DotArraySpecs
from .._lib import distributions as distr


class MainWidget(QWidget):

    def __init__(self, parent, settings, number, da_specs):
        super(MainWidget, self).__init__(parent)
        self.settings = settings

        self.initUI(number, da_specs)

    def initUI(self, number, da_specs):
        assert isinstance(da_specs, DotArraySpecs)
        assert isinstance(da_specs.diameter_distr, (distr.Beta, distr.Normal))

        self.btn_generate = QPushButton("Generate new array")

        self.number = LabeledNumberInput("Number", number)
        self.number2 = LabeledNumberInput("Number 2", 5)
        self.target_array_radius = LabeledNumberInput("Max radius", da_specs.target_array_radius)
        self.item_diameter_mean = LabeledNumberInput("Mean diameter", da_specs.diameter_distr.mu)
        self.item_diameter_std = LabeledNumberInput("Diameter range std", da_specs.diameter_distr.sigma)
        self.item_diameter_range = LabeledNumberInputTwoValues("Diameter range from",
                                                               value1=da_specs.diameter_distr.min_max[0],
                                                               value2=da_specs.diameter_distr.min_max[1])

        self.minimum_gap = LabeledNumberInput("Minimum gap", da_specs.minimum_gap)

        self.dot_colour = LabeledInput("Colour", text=self.settings.default_object_colour,
                                       case_sensitive=False)
        self.dot_colour2 = LabeledInput("Colour 2", text="expyriment_orange", case_sensitive=False)

        ctrl = QVBoxLayout()
        ctrl.addWidget(self.btn_generate)
        ctrl.addSpacing(10)
        ctrl.addLayout(self.number.layout())
        ctrl.addLayout(self.number2.layout())

        ctrl.addWidget(heading("Dot"))
        ctrl.addLayout(self.item_diameter_mean.layout())
        ctrl.addLayout(self.item_diameter_range.layout())
        ctrl.addLayout(self.item_diameter_std.layout())
        ctrl.addLayout(self.dot_colour.layout())
        ctrl.addLayout(self.dot_colour2.layout())

        ctrl.addWidget(heading("Array"))
        ctrl.addLayout(self.target_array_radius.layout())
        ctrl.addLayout(self.minimum_gap.layout())
        ctrl.addStretch(1)

        # Add text field
        fields = QVBoxLayout()
        self._text_field = QTextBrowser(self)
        self._text_field.setFont(QFont("Courier New", 8))
        self._text_field.setReadOnly(True)
        self._text_field.setWordWrapMode(False)
        self._text_field.setVerticalScrollBar(QScrollBar())
        self._text_field.setHorizontalScrollBar(QScrollBar())
        self.picture_field = QLabel(self)
        self.resize_fields(400)
        fields.addWidget(self.picture_field)
        fields.addWidget(self._text_field)

        hlayout = QHBoxLayout()
        hlayout.addLayout(ctrl)
        hlayout.addLayout(fields)
        self.setLayout(hlayout)

        self.updateUI()

    def updateUI(self):
        self.dot_colour2.setVisible(self.settings.bicoloured.isChecked())
        self.number2.setVisible(self.settings.bicoloured.isChecked())


    def resize_fields(self, width, text_height=150, minium_text_width=300):
        self.picture_field.setFixedSize(width, width)
        if width < minium_text_width:
            width = minium_text_width
        self._text_field.setMinimumSize(width, text_height)

    def text_clear(self):
        self._text_field.clear()

    def text_out(self, txt):
        self._text_field.append(txt)

    def text_error_feedback(self, error):
        self.text_clear()
        self.text_out("ERROR\n=====\n")
        txt = '<span style=\" color: #ff0000;\">%s</span>' % str(error)
        self.text_out(txt)
        self.updateUI()
