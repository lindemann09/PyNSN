"""
"""
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QCheckBox, \
    QDialogButtonBox, QComboBox, QHBoxLayout

from . import misc
from .._lib import constants
from .._visual_properties import flags, fit



class AdaptPropertyDialog(QDialog):

    def __init__(self, parent, properties):
        super(AdaptPropertyDialog, self).__init__(parent)

        self.setWindowTitle("Adapt Dot Array Property")
        self.properties = properties
        self._selection = None
        self.comboBox = QComboBox(self)
        for feat in flags:
            self.comboBox.addItem(feat.label())

        self.comboBox.activated[str].connect(self.choice)

        self._num_input = misc.NumberInput(width_edit=150, value=0)
        self.choice(flags.AV_DOT_DIAMETER)



        vlayout = QVBoxLayout(self)
        hlayout = QHBoxLayout()

        hlayout.addWidget(self.comboBox)
        hlayout.addWidget(self._num_input.edit)

        # OK and Cancel buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vlayout.addLayout(hlayout)
        vlayout.addWidget(buttons)

    def choice(self, selection):
        for feat in flags:
            if selection == feat:
                self._num_input.value = self.properties[feat.label()]
                self._selection = feat


    def current_state(self):
        return self._selection, self._num_input.value

    @staticmethod
    def get_response(parent, prop):
        """return the property to be adapted"""

        dialog = AdaptPropertyDialog(parent, prop)
        result = dialog.exec_()
        if result == QDialog.Accepted:
            return dialog.current_state()
        else:
            return None


class SettingsDialog(QDialog):

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
        self.colour_convex_hull_positions = misc.LabeledInput("Colour field position",
                                                              text=image_colours.field_area_positions.colour,
                                                              case_sensitive=False)
        self.colour_convex_hull_dots = misc.LabeledInput("Colour field area",
                                                         text=image_colours.field_area.colour,
                                                         case_sensitive=False)
        self.antialiasing = QCheckBox("Antialiasing")
        self.antialiasing.setChecked(True)

        self.bicoloured = QCheckBox("bicoloured")
        self.bicoloured.setChecked(False)

        self.default_object_colour = image_colours.default_object_colour

        vlayout = QVBoxLayout()
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
        buttons = QDialogButtonBox(QDialogButtonBox.Ok, Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)

        vlayout.addWidget(buttons)
        self.setLayout(vlayout)


class SequenceDialog(QDialog):
    extra_space = 50
    sequence_range = [10, 100]
    spacing_precision = constants.DEFAULT_FIT_SPACING_PRECISION
    adapt_FA2TA_ratio = constants.DEFAULT_FIT_FA2TA_RATIO


    def __init__(self, parent):

        super(SequenceDialog, self).__init__(parent)

        self.adapt_methods = []

        self.setWindowTitle("Sequence Dialog")

        self.adapt_diameter = QCheckBox(flags.AV_DOT_DIAMETER.label())
        self.adapt_av_perimeter= QCheckBox(flags.AV_PERIMETER.label())
        self.adapt_av_area = QCheckBox(flags.AV_SURFACE_AREA.label())
        self.adapt_area = QCheckBox(flags.TOTAL_SURFACE_AREA.label())
        self.adapt_total_perimeter = QCheckBox(flags.TOTAL_PERIMETER.label())
        self.adapt_coverage = QCheckBox(flags.COVERAGE.label())
        self.adapt_sparsity = QCheckBox(flags.SPARSITY.label())

        self.adapt_convex_hull = QCheckBox(flags.FIELD_AREA.label())
        self.adapt_size = QCheckBox(flags.LOG_SIZE.label())
        self.adapt_spacing = QCheckBox(flags.LOG_SPACING.label())
        self.adapt_spacing_presision = misc.LabeledNumberInput("Convex_hull presision",
                                                               value=SequenceDialog.spacing_precision,
                                                               integer_only=False)
        self.adapt_fa2ta = misc.LabeledNumberInput("Ratio convex_hull/area",
                                                   value=SequenceDialog.adapt_FA2TA_ratio,
                                                   integer_only=False, min=0, max=1)
        self.adapt_range = misc.LabeledNumberInputTwoValues("Sequence Range",
                                                            value1=SequenceDialog.sequence_range[0],
                                                            value2=SequenceDialog.sequence_range[1])
        self.adapt_extra_space = misc.LabeledNumberInput("Extra space",
                                                         value=SequenceDialog.extra_space,
                                                         integer_only=True, min=0)

        self.adapt_area.toggled.connect(self.ui_update)
        self.adapt_convex_hull.toggled.connect(self.ui_update)
        self.adapt_diameter.toggled.connect(self.ui_update)
        self.adapt_av_area.toggled.connect(self.ui_update)
        self.adapt_av_perimeter.toggled.connect(self.ui_update)
        self.adapt_total_perimeter.toggled.connect(self.ui_update)
        self.adapt_size.toggled.connect(self.ui_update)
        self.adapt_spacing.toggled.connect(self.ui_update)
        self.adapt_coverage.toggled.connect(self.ui_update)
        self.adapt_sparsity.toggled.connect(self.ui_update)
        self.adapt_spacing_presision.edit.editingFinished.connect(self.ui_update)
        self.adapt_fa2ta.edit.editingFinished.connect(self.ui_update)

        # OK and Cancel buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        vlayout = QVBoxLayout()
        vlayout.addWidget(misc.heading("Adapting"))
        vlayout.addWidget(self.adapt_diameter)
        vlayout.addWidget(self.adapt_av_perimeter)
        vlayout.addWidget(self.adapt_av_area)
        vlayout.addWidget(self.adapt_area)
        vlayout.addWidget(self.adapt_total_perimeter)
        vlayout.addWidget(self.adapt_convex_hull)
        vlayout.addWidget(self.adapt_coverage)
        vlayout.addWidget(self.adapt_sparsity)
        vlayout.addSpacing(10)
        vlayout.addWidget(self.adapt_size)
        vlayout.addWidget(self.adapt_spacing)
        vlayout.addSpacing(10)
        vlayout.addWidget(misc.heading("Adapting parameter"))
        vlayout.addLayout(self.adapt_spacing_presision.layout())
        vlayout.addLayout(self.adapt_fa2ta.layout())
        vlayout.addLayout(self.adapt_range.layout())
        vlayout.addSpacing(10)
        vlayout.addLayout(self.adapt_extra_space.layout())
        vlayout.addSpacing(10)

        vlayout.addWidget(buttons)
        self.setLayout(vlayout)
        self.ui_update()

    def ui_update(self):
        # get methods
        selected = []
        all = [flags.AV_DOT_DIAMETER]
        if self.adapt_diameter.isChecked():
            selected.append(all[-1])

        all.append(flags.AV_SURFACE_AREA)
        if self.adapt_av_area.isChecked():
            selected.append(all[-1])

        all.append(flags.AV_PERIMETER)
        if self.adapt_av_perimeter.isChecked():
            selected.append(all[-1])

        all.append(flags.TOTAL_PERIMETER)
        if self.adapt_total_perimeter.isChecked():
            selected.append(all[-1])

        all.append(flags.TOTAL_SURFACE_AREA)
        if self.adapt_area.isChecked():
            selected.append(all[-1])

        all.append(flags.FIELD_AREA)
        if self.adapt_convex_hull.isChecked():
            selected.append(all[-1])

        all.append(flags.COVERAGE)
        if self.adapt_coverage.isChecked():
            selected.append(all[-1])

        all.append(flags.SPARSITY)
        if self.adapt_sparsity.isChecked():
            selected.append(all[-1])

        all.append(flags.LOG_SIZE)
        if self.adapt_size.isChecked():
            selected.append(all[-1])

        all.append(flags.LOG_SPACING)
        if self.adapt_spacing.isChecked():
            selected.append(all[-1])

        self.adapt_diameter.setEnabled(True)
        self.adapt_av_perimeter.setEnabled(True)
        self.adapt_av_area.setEnabled(True)
        self.adapt_area.setEnabled(True)
        self.adapt_total_perimeter.setEnabled(True)
        self.adapt_coverage.setEnabled(True)
        self.adapt_sparsity.setEnabled(True)
        self.adapt_convex_hull.setEnabled(True)
        self.adapt_spacing.setEnabled(True)
        self.adapt_size.setEnabled(True)

        for x in all:
            if x not in selected:
                # test dependency of non-selected item, x, from any selected
                check = [s.is_dependent_from(x) for s in selected]
                if sum(check) > 0:  # any dependency
                    if x == flags.AV_DOT_DIAMETER:
                        self.adapt_diameter.setEnabled(False)
                        self.adapt_diameter.setChecked(False)
                    elif x == flags.AV_PERIMETER:
                        self.adapt_av_perimeter.setEnabled(False)
                        self.adapt_av_perimeter.setChecked(False)
                    elif x == flags.AV_SURFACE_AREA:
                        self.adapt_av_area.setEnabled(False)
                        self.adapt_av_area.setChecked(False)
                    elif x == flags.TOTAL_SURFACE_AREA:
                        self.adapt_area.setEnabled(False)
                        self.adapt_area.setChecked(False)
                    elif x == flags.TOTAL_PERIMETER:
                        self.adapt_total_perimeter.setEnabled(False)
                        self.adapt_total_perimeter.setChecked(False)
                    elif x == flags.FIELD_AREA:
                        self.adapt_convex_hull.setEnabled(False)
                        self.adapt_convex_hull.setChecked(False)
                    elif x == flags.COVERAGE:
                        self.adapt_coverage.setEnabled(False)
                        self.adapt_coverage.setChecked(False)
                    elif x == flags.SPARSITY:
                        self.adapt_sparsity.setEnabled(False)
                        self.adapt_sparsity.setChecked(False)
                    elif x == flags.LOG_SIZE:
                        self.adapt_size.setEnabled(False)
                        self.adapt_size.setChecked(False)
                    elif x == flags.LOG_SPACING:
                        self.adapt_spacing.setEnabled(False)
                        self.adapt_spacing.setChecked(False)

        self.adapt_methods = selected

    @staticmethod
    def get_response(parent):

        dialog = SequenceDialog(parent)
        result = dialog.exec_()
        if result == QDialog.Accepted:
            SequenceDialog.extra_space = dialog.adapt_extra_space.value
            SequenceDialog.spacing_precision = dialog.adapt_spacing_presision.value
            SequenceDialog.adapt_FA2TA_ratio = dialog.adapt_fa2ta.value
            SequenceDialog.sequence_range = [dialog.adapt_range.value1, dialog.adapt_range.value2]
            return (dialog.adapt_methods,
                    SequenceDialog.sequence_range,
                    SequenceDialog.extra_space)
        else:
            return (None, None, None)
