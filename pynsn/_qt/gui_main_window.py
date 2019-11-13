
import json
import os
from PyQt4 import QtGui
from PIL.ImageQt import ImageQt
from .._lib import random_dot_array
from .._lib import dot_array_sequence
from .._lib._logging import LogFile
from .._lib import visual_features as vf
from .._lib import _colour
from .._lib._item_attributes import ItemAttributes
from .._lib import pil_image
from .main_widget import MainWidget
from . import dialogs
from .sequence_display import SequenceDisplay


DEFAULT_ARRAY = (40, random_dot_array.Specs(target_area_radius=200,
                                            item_colour="lime",
                                            item_diameter_mean=15,
                                            item_diameter_range=[5, 40],
                                            item_diameter_std=8,
                                            minimum_gap=2),
                 _colour .ImageColours(target_area="#3e3e3e",
                                        field_area=None,
                                        field_area_outer=None,
                                        center_of_mass=None,
                                        center_of_outer_positions=None,
                                        background="gray"))

ICON = (11, random_dot_array.Specs(target_area_radius=200,
                                   item_colour="lime",
                                   item_diameter_mean=35,
                                   item_diameter_range=[5, 80],
                                   item_diameter_std=20),
        _colour .ImageColours(target_area="#3e3e3e",
                               field_area=None,
                               field_area_outer="expyriment_orange",
                               center_of_mass=None,
                               center_of_outer_positions=None,
                               background=None))


class GUIMainWindow(QtGui.QMainWindow):

    def __init__(self):

        super(GUIMainWindow, self).__init__()

        self._image = None
        self.dot_array = None
        self.set_loging(False)  # todo checkbox in settings
        self.settings = dialogs.SettingsDialog(self, image_colours=DEFAULT_ARRAY[2])
        self.initUI()
        self.show()

    def set_loging(self, onoff):

        if onoff:
            self.logger = LogFile(log_filename="log/_qt",
                                  override_log_files=True,
                                  log_colours=False)
        else:
            self.logger = None

    def initUI(self):

        # menues
        exitAction = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)

        settingsAction = QtGui.QAction('&Settings', self)
        settingsAction.triggered.connect(self.action_settings)

        saveAction = QtGui.QAction('&Save current stimulus', self)
        saveAction.setShortcut('Ctrl+S')
        saveAction.triggered.connect(self.save_array)


        printxyAction = QtGui.QAction('&Print array', self)
        printxyAction.triggered.connect(self.action_print_xy)
        printparaAction = QtGui.QAction('&Print parameter', self)
        printparaAction.triggered.connect(self.action_print_para)

        matchAction = QtGui.QAction('&Match property', self)
        matchAction.triggered.connect(self.action_match)

        sequenceAction = QtGui.QAction('&Make sequence', self)
        sequenceAction.triggered.connect(self.action_make_sequence)

        aboutAction = QtGui.QAction('&About', self)

        # self.statusBar()
        menubar = self.menuBar()

        arrayMenu = menubar.addMenu('&Array')
        arrayMenu.addAction(settingsAction)
        arrayMenu.addAction(saveAction)
        arrayMenu.addSeparator()
        arrayMenu.addAction(exitAction)

        toolMenu = menubar.addMenu('&Tools')
        toolMenu.addAction(sequenceAction)
        toolMenu.addAction(matchAction)
        toolMenu.addSeparator()
        toolMenu.addAction(printparaAction)
        toolMenu.addAction(printxyAction)

        aboutMenu = menubar.addMenu('&About')
        aboutMenu.addAction(aboutAction)



        # main widget
        self.main_widget = MainWidget(self, self.settings, DEFAULT_ARRAY[0], DEFAULT_ARRAY[1])
        self.setCentralWidget(self.main_widget)
        self.main_widget.btn_generate.clicked.connect(self.action_generate_btn)
        self.main_widget.dot_colour.edit.editingFinished.connect(self.action_dot_colour_change)
        self.main_widget.slider.sliderReleased.connect(self.action_slider_released)

        self.move(300, -300)
        self.setWindowTitle('PyNSN GUI')

        # ICON
        colours = ICON[2]
        self._image = pil_image.create(
                        dot_array=random_dot_array.create(n_dots=ICON[0], specs = ICON[1]),
                        colours=colours, antialiasing=True)

        self.setWindowIcon(QtGui.QIcon(self.pixmap()))
        self._image = None

        self.action_generate_btn()

    def make_new_array(self):

        try:
            self.dot_array = random_dot_array.create(n_dots=self.get_number(),
                                                     specs=self.get_specs(),
                                                     logger=self.logger)
        except (RuntimeError, ValueError) as error:
            self.main_widget.text_error_feedback(error)
            raise error

        if self.settings.bicoloured.isChecked():
            data_array2 = random_dot_array.create(n_dots=self.main_widget.number2.value,
                                                  specs=self.get_specs(),
                                                  occupied_space=self.dot_array,
                                                  logger=self.logger)
            data_array2.set_attributes(
                ItemAttributes(colour=self.main_widget.dot_colour2.text))
            self.dot_array.join(data_array2, realign=False)

        self.dot_array.round(decimals=self.settings.rounding_decimals.value)
        self._image = None

    def image(self):
        if self._image is not None:
            return self._image
        else:
            para = self.get_image_colours()
            image_colours = _colour.ImageColours(
                target_area=para.target_area,
                field_area=para.field_area,
                field_area_outer=para.field_area_outer,
                center_of_mass=para.center_of_mass,
                center_of_outer_positions=para.center_of_outer_positions,
                background=para.background)

            self._image = pil_image.create(dot_array=self.dot_array,
                                           colours=image_colours,
                                           antialiasing=self.settings.antialiasing.isChecked())
                                            # todo maybe: gabor_filter=ImageFilter.GaussianBlur

            return self._image

    def get_number(self):
        return self.main_widget.number.value

    def get_specs(self):

        try:
            colour_dot = _colour.Colour(self.main_widget.dot_colour.text)
        except:
            colour_dot = DEFAULT_ARRAY[1].item_attributes.colour
            self.main_widget.dot_colour.text = colour_dot

        return random_dot_array.Specs(target_area_radius=self.main_widget.target_array_radius.value,
                                      item_colour=colour_dot,
                                      item_diameter_mean=self.main_widget.item_diameter_mean.value,
                                      item_diameter_range=[self.main_widget.item_diameter_range.value1,
                                                      self.main_widget.item_diameter_range.value2],
                                      item_diameter_std=self.main_widget.item_diameter_std.value,
                                      minimum_gap=self.main_widget.minimum_gap.value)

    def get_image_colours(self):
        # check colour input

        try:
            colour_area = _colour.Colour(self.settings.colour_area.text)
        except:
            colour_area = None
            self.settings.colour_area.text = "None"
        try:
            colour_convex_hull_positions = _colour.Colour(
                self.settings.colour_convex_hull_positions.text)
        except:
            colour_convex_hull_positions = None
            self.settings.colour_convex_hull_positions.text = "None"
        try:
            colour_convex_hull_dots = _colour.Colour(
                self.settings.colour_convex_hull_dots.text)
        except:
            colour_convex_hull_dots = None
            self.settings.colour_convex_hull_dots.text = "None"
        try:
            colour_background = _colour.Colour(
                self.settings.colour_background.text)
        except:
            colour_background = None
            self.settings.colour_background.text = "None"

        return _colour.ImageColours(target_area=colour_area,
                                      field_area=colour_convex_hull_positions,
                                      field_area_outer=colour_convex_hull_dots,
                                      center_of_mass=None,
                                      center_of_outer_positions=None,
                                      background=colour_background)

    def pixmap(self):
        return QtGui.QPixmap.fromImage(ImageQt(self.image()))

    def show_current_image(self, remake_image=False):
        """"""
        if remake_image:
            self._image = None
        w = self.get_specs().target_array_radius * 2
        self.main_widget.resize_fields(width=w, text_height=150)
        self.main_widget.picture_field.setPixmap(self.pixmap())
        self.main_widget.adjustSize()
        self.adjustSize()

    def action_generate_btn(self):
        """"""
        self.make_new_array()
        self.show_current_image(remake_image=True)
        self.write_properties()
        self.main_widget.updateUI()

    def write_properties(self, clear_field=True):
        txt = self.dot_array.feature.get_features_text(extended_format=True,
                                                       with_hash=True)
        if self.settings.bicoloured.isChecked():
            for da in self.dot_array.split_array_by_colour():
                txt += "Colour {}\n".format(da.get_colours()[0])
                txt += da.feature.get_features_text(extended_format=True,
                                            with_hash=False)
        if clear_field:
            self.main_widget.text_clear()

        self.main_widget.text_out(txt)

    def action_print_xy(self):
        """"""
        txt = self.dot_array.get_csv(hash_column=False, num_idx_column=False, colour_column=True)
        self.main_widget.text_out(txt)

    def action_print_para(self):
        d = {'number': self.get_number()}
        d['image_parameter'] = self.get_image_colours().as_dict()
        d['randomization'] = self.get_specs().as_dict()
        self.main_widget.text_out("# parameter\n" + json.dumps(d, indent=2))

    def save_array(self):
        """"""

        filename, extension = QtGui.QFileDialog.getSaveFileNameAndFilter(
            self,
            'Save file',
            "",
            "Image PNG File (*.png);; Image BMP File (*.bmp);; JSON File (.json)")

        if len(filename)>0:
            filename = os.path.abspath(filename)

            ext = extension.split("(")[1].replace(")", "")\
                                         .replace("*", "").strip()
            if not filename.endswith(ext):
                filename = filename + ext

            if ext == ".json":
                self.dot_array.save(filename)

            elif ext == ".png" or ext == ".bmp":
                self.image().save(filename)


    def action_match(self):
        """"""
        prop = self.dot_array.feature.get_features_dict()
        response = dialogs.MatchPropertyDialog.get_response(self, prop)  #
        if response is not None:
            self.dot_array.match(response)
            if isinstance(response, vf.SIZE_FEATURES):
                self.dot_array.center_array()
                self.dot_array.realign()
            self.show_current_image(remake_image=True)
            self.write_properties()
            self.main_widget.updateUI()


    def action_settings(self):
        """"""
        result = self.settings.exec_()
        self.main_widget.updateUI()
        self.show_current_image(remake_image=True)

    def action_dot_colour_change(self):
        """"""
        self.dot_array.set_attributes(self.get_specs().item_attributes)
        self.show_current_image(remake_image=True)

    def action_slider_released(self):
        """"""
        change = self.main_widget.number.value - self.dot_array.feature.numerosity
        self.dot_array = self.dot_array.number_deviant(change)
        self.show_current_image(remake_image=True)
        self.write_properties()
        # todo slider does not work correctly for multi colour arrays

    def action_make_sequence(self):
        match_methods, match_range, extra_space = dialogs.SequenceDialog.get_response(self)

        d = {"match range": match_range, "extra_space": extra_space}
        d["match_methods"] = list(map(lambda x: x.get_features_dict(), match_methods))  # TODO <-- check ERROR
        self.main_widget.text_out("# Sequence\n" + \
                                           json.dumps(d))

        if match_methods is not None:
            # print("processing")
            sequence = dot_array_sequence.create(
                reference_dot_array=self.dot_array,
                              match_properties=match_methods,
                              min_max_numerosity=match_range,
                              extra_space=extra_space,
                              logger=self.logger)
            SequenceDisplay(self, da_sequence=sequence,
                            start_numerosity=self.dot_array.feature.numerosity,
                            image_colours=self.get_image_colours(),
                            antialiasing=self.settings.antialiasing.isChecked()).exec_()


