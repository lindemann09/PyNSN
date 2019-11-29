"""
"""
from PyQt4 import QtGui, QtCore
from . import misc
from .._lib._visual_features import Features
from .._lib._match import DEFAULT_SPACING_PRECISION, DEFAULT_MATCH_FA2TA_RATIO

class MatchPropertyDialog(QtGui.QDialog):

    def __init__(self, parent, properties):
        super(MatchPropertyDialog, self).__init__(parent)

        self.setWindowTitle("Match Dot Array Property")
        self.features = properties
        self._selection = None
        self.comboBox = QtGui.QComboBox(self)
        for feat in Features.ALL_FEATURES:
            self.comboBox.addItem(feat)

        self.comboBox.activated[str].connect(self.choice)

        self._num_input = misc.NumberInput(width_edit=150, value=0)
        self.choice(Features.ITEM_DIAMETER)



        vlayout = QtGui.QVBoxLayout(self)
        hlayout = QtGui.QHBoxLayout()

        hlayout.addWidget(self.comboBox)
        hlayout.addWidget(self._num_input.edit)

        # OK and Cancel buttons
        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vlayout.addLayout(hlayout)
        vlayout.addWidget(buttons)

    def choice(self, selection):

        for feat in Features.ALL_FEATURES:
            if selection == feat:
                self._num_input.value = self.features[feat]
                self._selection = feat


    def current_state(self):
        return self._selection, self._num_input.value

    @staticmethod
    def get_response(parent, prop):
        """return the feature to be matched"""

        dialog = MatchPropertyDialog(parent, prop)
        result = dialog.exec_()
        if result == QtGui.QDialog.Accepted:
            return dialog.current_state()
        else:
            return None


class SettingsDialog(QtGui.QDialog):

    def __init__(self, parent, image_colours):
        super(SettingsDialog, self).__init__(parent)


        self.setWindowTitle("Dot Array Property")

        self.rounding_decimals = misc.LabeledNumberInput(
                            label="Rounding decimals",
                            value=0, integer_only = True,
                            min = 0, max = 8)

        self.colour_area = misc.LabeledInput("Target Area",
                                             text=image_colours.target_area.colour,
                                             case_sensitive=False)
        self.colour_background = misc.LabeledInput("Background",
                                                   text=image_colours.background.colour,
                                                   case_sensitive=False)
        self.colour_convex_hull_positions = misc.LabeledInput("Colour field area",
                                                              text=image_colours.field_area.colour,
                                                              case_sensitive=False)
        self.colour_convex_hull_dots = misc.LabeledInput("Colour field area outer",
                                                         text=image_colours.field_area_outer.colour,
                                                         case_sensitive=False)
        self.antialiasing = QtGui.QCheckBox("Antialiasing")
        self.antialiasing.setChecked(True)

        self.bicoloured = QtGui.QCheckBox("bicoloured")
        self.bicoloured.setChecked(False)



        vlayout = QtGui.QVBoxLayout()
        vlayout.addLayout(self.rounding_decimals.layout())

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
    extra_space = 50
    sequence_range = [10, 100]
    spacing_precision = DEFAULT_SPACING_PRECISION
    match_FA2TA_ratio = DEFAULT_MATCH_FA2TA_RATIO

    def __init__(self, parent):

        super(SequenceDialog, self).__init__(parent)

        self.match_methods = []

        self.setWindowTitle("Sequence Dialog")

        self.match_diameter = QtGui.QCheckBox(Features.ITEM_DIAMETER)
        self.match_item_perimeter= QtGui.QCheckBox(Features.ITEM_PERIMETER)
        self.match_item_area = QtGui.QCheckBox(Features.ITEM_SURFACE_AREA)
        self.match_area = QtGui.QCheckBox(Features.TOTAL_SURFACE_AREA)
        self.match_total_perimeter = QtGui.QCheckBox(Features.TOTAL_PERIMETER)
        self.match_coverage = QtGui.QCheckBox(Features.COVERAGE)
        self.match_sparsity = QtGui.QCheckBox(Features.SPARSITY)

        self.match_convex_hull = QtGui.QCheckBox(Features.FIELD_AREA)
        self.match_size = QtGui.QCheckBox(Features.LOG_SIZE)
        self.match_spacing = QtGui.QCheckBox(Features.LOG_SPACING)
        self.match_spacing_presision = misc.LabeledNumberInput("Convex_hull presision",
                                                               value=SequenceDialog.spacing_precision,
                                                               integer_only=False)
        self.match_fa2ta = misc.LabeledNumberInput("Ratio convex_hull/area",
                                                   value=SequenceDialog.match_FA2TA_ratio,
                                                   integer_only=False, min=0, max=1)
        self.match_range = misc.LabeledNumberInputTwoValues("Sequence Range",
                                                            value1=SequenceDialog.sequence_range[0],
                                                            value2=SequenceDialog.sequence_range[1])
        self.match_extra_space = misc.LabeledNumberInput("Extra space",
                                                         value=SequenceDialog.extra_space,
                                                         integer_only=True, min=0)

        self.match_area.toggled.connect(self.ui_update)
        self.match_convex_hull.toggled.connect(self.ui_update)
        self.match_diameter.toggled.connect(self.ui_update)
        self.match_item_area.toggled.connect(self.ui_update)
        self.match_item_perimeter.toggled.connect(self.ui_update)
        self.match_total_perimeter.toggled.connect(self.ui_update)
        self.match_size.toggled.connect(self.ui_update)
        self.match_spacing.toggled.connect(self.ui_update)
        self.match_coverage.toggled.connect(self.ui_update)
        self.match_sparsity.toggled.connect(self.ui_update)
        self.match_spacing_presision.edit.editingFinished.connect(self.ui_update)
        self.match_fa2ta.edit.editingFinished.connect(self.ui_update)

        # OK and Cancel buttons
        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vlayout = QtGui.QVBoxLayout()
        vlayout.addWidget(misc.heading("Matching"))
        vlayout.addWidget(self.match_diameter)
        vlayout.addWidget(self.match_item_perimeter)
        vlayout.addWidget(self.match_item_area)
        vlayout.addWidget(self.match_area)
        vlayout.addWidget(self.match_total_perimeter)
        vlayout.addWidget(self.match_convex_hull)
        vlayout.addWidget(self.match_coverage)
        vlayout.addWidget(self.match_sparsity)
        vlayout.addSpacing(10)
        vlayout.addWidget(self.match_size)
        vlayout.addWidget(self.match_spacing)
        vlayout.addSpacing(10)
        vlayout.addWidget(misc.heading("Matching parameter"))
        vlayout.addLayout(self.match_spacing_presision.layout())
        vlayout.addLayout(self.match_fa2ta.layout())
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
        all = [Features.ITEM_DIAMETER]
        if self.match_diameter.isChecked():
            selected.append(all[-1])

        all.append(Features.ITEM_SURFACE_AREA)
        if self.match_item_area.isChecked():
            selected.append(all[-1])

        all.append(Features.ITEM_PERIMETER)
        if self.match_item_perimeter.isChecked():
            selected.append(all[-1])

        all.append(Features.TOTAL_PERIMETER)
        if self.match_total_perimeter.isChecked():
            selected.append(all[-1])

        all.append(Features.TOTAL_SURFACE_AREA)
        if self.match_area.isChecked():
            selected.append(all[-1])

        all.append(Features.FIELD_AREA)
        if self.match_convex_hull.isChecked():
            selected.append(all[-1])

        all.append(Features.COVERAGE)
        if self.match_coverage.isChecked():
            selected.append(all[-1])

        all.append(Features.SPARSITY)
        if self.match_sparsity.isChecked():
            selected.append(all[-1])

        all.append(Features.LOG_SIZE)
        if self.match_size.isChecked():
            selected.append(all[-1])

        all.append(Features.LOG_SPACING)
        if self.match_spacing.isChecked():
            selected.append(all[-1])

        self.match_diameter.setEnabled(True)
        self.match_item_perimeter.setEnabled(True)
        self.match_item_area.setEnabled(True)
        self.match_area.setEnabled(True)
        self.match_total_perimeter.setEnabled(True)
        self.match_coverage.setEnabled(True)
        self.match_sparsity.setEnabled(True)
        self.match_convex_hull.setEnabled(True)
        self.match_spacing.setEnabled(True)
        self.match_size.setEnabled(True)

        for x in all:
            if x not in selected:
                # test dependency of non-selected item, x, from any selected
                check = map(lambda s: Features.are_dependent(x, s), selected)
                if sum(check) > 0:  # any dependency
                    if x == Features.ITEM_DIAMETER:
                        self.match_diameter.setEnabled(False)
                        self.match_diameter.setChecked(False)
                    elif x == Features.ITEM_PERIMETER:
                        self.match_item_perimeter.setEnabled(False)
                        self.match_item_perimeter.setChecked(False)
                    elif x == Features.ITEM_SURFACE_AREA:
                        self.match_item_area.setEnabled(False)
                        self.match_item_area.setChecked(False)
                    elif x == Features.TOTAL_SURFACE_AREA:
                        self.match_area.setEnabled(False)
                        self.match_area.setChecked(False)
                    elif x == Features.TOTAL_PERIMETER:
                        self.match_total_perimeter.setEnabled(False)
                        self.match_total_perimeter.setChecked(False)
                    elif x == Features.FIELD_AREA:
                        self.match_convex_hull.setEnabled(False)
                        self.match_convex_hull.setChecked(False)
                    elif x == Features.COVERAGE:
                        self.match_coverage.setEnabled(False)
                        self.match_coverage.setChecked(False)
                    elif x == Features.SPARSITY:
                        self.match_sparsity.setEnabled(False)
                        self.match_sparsity.setChecked(False)
                    elif x == Features.LOG_SIZE:
                        self.match_size.setEnabled(False)
                        self.match_size.setChecked(False)
                    elif x == Features.LOG_SPACING:
                        self.match_spacing.setEnabled(False)
                        self.match_spacing.setChecked(False)

        self.match_methods = selected

    @staticmethod
    def get_response(parent):

        dialog = SequenceDialog(parent)
        result = dialog.exec_()
        if result == QtGui.QDialog.Accepted:
            SequenceDialog.extra_space = dialog.match_extra_space.value
            SequenceDialog.spacing_precision = dialog.match_spacing_presision.value
            SequenceDialog.match_FA2TA_ratio = dialog.match_fa2ta.value
            SequenceDialog.sequence_range = [dialog.match_range.value1, dialog.match_range.value2]

            return (dialog.match_methods, SequenceDialog.sequence_range,
                    SequenceDialog.extra_space)
        else:
            return (None, None, None)
