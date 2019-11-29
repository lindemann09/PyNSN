"""
"""
from PyQt4 import QtGui, QtCore
from .misc import heading, LabeledNumberInput, LabeledNumberInputTwoValues, LabeledInput

class MainWidget(QtGui.QWidget):

    def __init__(self, parent, settings, number, generator):
        super(MainWidget, self).__init__(parent)
        self.settings = settings
        self.initUI(number, generator)

    def initUI(self, number, generator):
        self.btn_generate = QtGui.QPushButton("Generate new array")

        self.number = LabeledNumberInput("Number", number)
        self.number2 = LabeledNumberInput("Number 2", 0)

        self.target_array_radius = LabeledNumberInput("Max radius", generator.target_array_radius)
        self.item_diameter_mean = LabeledNumberInput("Mean diameter", generator.item_diameter_mean)
        self.item_diameter_std = LabeledNumberInput("Diameter range std", generator.item_diameter_std)
        self.item_diameter_range = LabeledNumberInputTwoValues("Diameter range from",
                                                               value1=generator.item_diameter_range[0],
                                                               value2=generator.item_diameter_range[1])

        self.minimum_gap = LabeledNumberInput("Minimum gap", generator.minimum_gap)

        self.dot_colour = LabeledInput("Colour", text=generator.item_attributes.colour, case_sensitive=False)
        self.dot_colour2 = LabeledInput("Colour 2", text="skyblue", case_sensitive=False)

        ctrl = QtGui.QVBoxLayout()
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
        fields = QtGui.QVBoxLayout()
        self._text_field = QtGui.QTextBrowser(self)
        self._text_field.setFont(QtGui.QFont("Courier New", 8))
        self._text_field.setReadOnly(True)
        self._text_field.setWordWrapMode(False)
        self._text_field.setVerticalScrollBar(QtGui.QScrollBar())
        self._text_field.setHorizontalScrollBar(QtGui.QScrollBar())
        self.picture_field = QtGui.QLabel(self)
        self.resize_fields(400)
        fields.addWidget(self.picture_field)
        fields.addWidget(self._text_field)

        hlayout = QtGui.QHBoxLayout()
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

