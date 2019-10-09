"""
"""
from PyQt4 import QtGui, QtCore
from . import misc
from pynsn._lib import visual_features as vf


class MatchPropertyDialog(QtGui.QDialog):

    def __init__(self, parent, properties):
        super(MatchPropertyDialog, self).__init__(parent)

        self.setWindowTitle("Match Dot Array Property")
        self.features = properties
        self._selection = None
        self.comboBox = QtGui.QComboBox(self)
        for feat in vf.ALL_VISUAL_FEATURES:
            self.comboBox.addItem(feat.label)

        self.comboBox.activated[str].connect(self.choice)

        self._num_input = misc.NumberInput(width_edit=150, value=0)
        self.choice(vf.ItemDiameter.label)



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

        for feat in vf.ALL_VISUAL_FEATURES:
            if selection == feat.label:
                self._num_input.value = self.features[feat.label]
                self._selection = feat()


    def current_state(self):
        self._selection.value = self._num_input.value
        return self._selection

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
    def __init__(self, parent, image_parameter):
        super(SettingsDialog, self).__init__(parent)

        self.setWindowTitle("Dot Array Property")

        self.colour_area = misc.LabeledInput("Traget Area",
                                             text=image_parameter.colour_target_area.colour,
                                             case_sensitive=False)
        self.colour_background = misc.LabeledInput("Background",
                                                   text=image_parameter.colour_background.colour,
                                                   case_sensitive=False)
        self.colour_convex_hull_positions = misc.LabeledInput("Colour field area",
                                                              text=image_parameter.colour_field_area.colour,
                                                              case_sensitive=False)
        self.colour_convex_hull_dots = misc.LabeledInput("Colour field area outer",
                                                         text=image_parameter.colour_field_area_outer.colour,
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
    extra_space = 50
    spacing_precision = vf.FieldArea().spacing_precision
    fieldarea2totalarea = vf.Coverage().match_ratio_fieldarea2totalarea
    sequence_range = [10, 100]

    def __init__(self, parent):

        super(SequenceDialog, self).__init__(parent)

        self.match_methods = []

        self.setWindowTitle("Sequence Dialog")

        self.match_diameter = QtGui.QCheckBox(vf.ItemDiameter.label)
        self.match_item_perimeter= QtGui.QCheckBox(vf.ItemPerimeter.label)
        self.match_item_area = QtGui.QCheckBox(vf.ItemSurfaceArea.label)
        self.match_area = QtGui.QCheckBox(vf.TotalSurfaceArea.label)
        self.match_total_perimeter = QtGui.QCheckBox(vf.TotalPerimeter.label)
        self.match_coverage = QtGui.QCheckBox(vf.Coverage.label)
        self.match_sparsity = QtGui.QCheckBox(vf.Sparsity.label)

        self.match_convex_hull = QtGui.QCheckBox(vf.FieldArea.label)
        self.match_size = QtGui.QCheckBox(vf.LogSize.label)
        self.match_spacing = QtGui.QCheckBox(vf.LogSpacing.label)
        self.match_spacing_presision = misc.LabeledNumberInput("Convex_hull presision",
                                                               value=SequenceDialog.spacing_precision,
                                                               integer_only=False)
        self.match_fa2ta = misc.LabeledNumberInput("Ratio convex_hull/area",
                                                   value=SequenceDialog.fieldarea2totalarea,
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
        all = [vf.ItemDiameter()]
        if self.match_diameter.isChecked():
            selected.append(all[-1])

        all.append(vf.ItemSurfaceArea())
        if self.match_item_area.isChecked():
            selected.append(all[-1])

        all.append(vf.ItemPerimeter())
        if self.match_item_perimeter.isChecked():
            selected.append(all[-1])

        all.append(vf.TotalPerimeter())
        if self.match_total_perimeter.isChecked():
            selected.append(all[-1])

        all.append(vf.TotalSurfaceArea())
        if self.match_area.isChecked():
            selected.append(all[-1])

        all.append(vf.FieldArea(spacing_precision=self.match_spacing_presision.value))
        if self.match_convex_hull.isChecked():
            selected.append(all[-1])

        all.append(vf.Coverage(match_ratio_fieldarea2totalarea=self.match_fa2ta.value,
                                 spacing_precision=self.match_spacing_presision.value))
        if self.match_coverage.isChecked():
            selected.append(all[-1])

        all.append(vf.Sparsity(spacing_precision=self.match_spacing_presision.value))
        if self.match_sparsity.isChecked():
            selected.append(all[-1])


        all.append(vf.LogSize())
        if self.match_size.isChecked():
            selected.append(all[-1])

        all.append(vf.LogSpacing(spacing_precision=self.match_spacing_presision.value))
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
                if sum(map(lambda s: s.is_dependent(x), selected)) > 0:  # any dependency
                    if isinstance(x, vf.ItemDiameter):
                        self.match_diameter.setEnabled(False)
                        self.match_diameter.setChecked(False)
                    if isinstance(x, vf.ItemPerimeter):
                        self.match_item_perimeter.setEnabled(False)
                        self.match_item_perimeter.setChecked(False)
                    if isinstance(x, vf.ItemSurfaceArea):
                        self.match_item_area.setEnabled(False)
                        self.match_item_area.setChecked(False)
                    if isinstance(x, vf.TotalSurfaceArea):
                        self.match_area.setEnabled(False)
                        self.match_area.setChecked(False)
                    if isinstance(x, vf.TotalPerimeter):
                        self.match_total_perimeter.setEnabled(False)
                        self.match_total_perimeter.setChecked(False)
                    if isinstance(x, vf.FieldArea):
                        self.match_convex_hull.setEnabled(False)
                        self.match_convex_hull.setChecked(False)
                    if isinstance(x, vf.Coverage):
                        self.match_coverage.setEnabled(False)
                        self.match_coverage.setChecked(False)
                    if isinstance(x, vf.Sparsity):
                        self.match_sparsity.setEnabled(False)
                        self.match_sparsity.setChecked(False)
                    if isinstance(x, vf.LogSize):
                        self.match_size.setEnabled(False)
                        self.match_size.setChecked(False)
                    if isinstance(x, vf.LogSpacing):
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
            SequenceDialog.fieldarea2totalarea = dialog.match_fa2ta.value
            SequenceDialog.sequence_range = [dialog.match_range.value1, dialog.match_range.value2]

            return (dialog.match_methods, SequenceDialog.sequence_range,
                    SequenceDialog.extra_space)
        else:
            return (None, None, None)
