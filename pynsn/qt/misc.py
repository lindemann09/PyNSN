"""
"""
from PyQt5.QtCore import Qt, QCoreApplication
from PyQt5.QtGui import QFont, QIntValidator, QDoubleValidator
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QProgressDialog


class LabeledInput(object):

    def __init__(self, label, text, width_label=180, width_edit=70, case_sensitive=True):

        self.label = QLabel(label)
        self.label.setFixedWidth(width_label)
        self.edit = QLineEdit()
        self.edit.setFixedWidth(width_edit)
        self.edit.setAlignment(Qt.AlignRight)
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
            layout = QVBoxLayout()
        else:
            layout = QHBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.edit)
        return layout

    def setVisible(self, bool):
        self.edit.setVisible(bool)
        self.label.setVisible(bool)


class LabeledNumberInput(LabeledInput):
    def __init__(self, label, value, width_label=180, width_edit=70,
                 integer_only=True, min=None, max=None):

        LabeledInput.__init__(self, label=label, text="", width_label=width_label,
                              width_edit=width_edit)
        self.integer_only = integer_only
        if integer_only:
            self.edit.setValidator(QIntValidator())
        else:
            self.edit.setValidator(QDoubleValidator())
        self.min = min
        self.max = max
        self.value = value

    @property
    def value(self):
        if self.integer_only:
            rtn = int(self.edit.text())
        else:
            rtn = float(self.edit.text())
        if self.min is not None and rtn < self.min:
            rtn = self.min
            self.value = rtn
        if self.max is not None and rtn > self.max:
            rtn = self.max
            self.value = rtn

        return rtn

    @value.setter
    def value(self, v):
        if self.min is not None and v < self.min:
            v = self.min
        if self.max is not None and v > self.max:
            v = self.max
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
    boldFont = QFont()
    boldFont.setBold(True)
    heading = QLabel(text)
    heading.setFont(boldFont)
    heading.setStyleSheet("QLabel { color : black; }")
    return heading


class NumberInput(object):

    def __init__(self, value, width_edit=70, integer_only=False, min=None, max=None):

        self.edit = QLineEdit()
        self.edit.setFixedWidth(width_edit)
        self.edit.setAlignment(Qt.AlignRight)
        self.integer_only = integer_only
        if integer_only:
            self.edit.setValidator(QIntValidator())
        else:
            self.edit.setValidator(QDoubleValidator())
        self.min = min
        self.max = max
        self.value = value

    @property
    def value(self):
        if self.integer_only:
            rtn = int(self.edit.text())
        else:
            rtn = float(self.edit.text())
        if self.min is not None and rtn < self.min:
            rtn = self.min
            self.value = rtn
        if self.max is not None and rtn > self.max:
            rtn = self.max
            self.value = rtn

        return rtn

    @value.setter
    def value(self, v):
        if self.min is not None and v < self.min:
            v = self.min
        if self.max is not None and v > self.max:
            v = self.max
        if self.integer_only:
            v = int(v)
        else:
            v = float(v)
        self.edit.setText(str(v))


def progressbar_iterator(iteratable, n_elements, label, win_title):
    # iterator function with prgress bar
    dialog = QProgressDialog()
    dialog.setMaximum(100)
    dialog.setLabelText(label)
    dialog.setWindowTitle(win_title)
    dialog.open()
    for cnt, item in enumerate(iteratable):
        p = (cnt * 100) // n_elements
        QCoreApplication.instance().processEvents()
        if dialog.wasCanceled():
            raise StopIteration
        dialog.setValue(p)
        yield (item)

    dialog.close()
