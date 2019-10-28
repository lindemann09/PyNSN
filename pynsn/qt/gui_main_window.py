
import yaml
import os
from PyQt4 import QtGui
from PIL.ImageQt import ImageQt
from ..lib.dot_array import DotArrayFactory
from ..lib.dot_array_sequence import generate_da_sequence
from ..lib.logging import LogFile
from ..lib.colour import Colour
from ..lib import pil_image
from .main_widget import MainWidget
from . import dialogs
from .sequence_display import SequenceDisplay

DEFAULT_ARRAY = (40, DotArrayFactory(target_area_radius=200,
                                     item_colour="lime",
                                     item_diameter_mean=15,
                                     item_diameter_range=[5, 40],
                                     item_diameter_std=8,
                                     minimum_gap=2),
                 pil_image.PILImagePlotter(colour_target_area="#3e3e3e",
                                           colour_field_area=None,
                                           colour_field_area_outer=None,
                                           colour_center_of_mass=None,
                                           colour_center_of_outer_positions=None,
                                           antialiasing=True,
                                           colour_background="gray"))

ICON = (11, DotArrayFactory(target_area_radius=200,
                            item_colour="lime",
                            item_diameter_mean=35,
                            item_diameter_range=[5, 80],
                            item_diameter_std=20),
        pil_image.PILImagePlotter(colour_target_area="#3e3e3e",
                                  colour_field_area=None,
                                  colour_field_area_outer="expyriment_orange",
                                  colour_center_of_mass=None,
                                  colour_center_of_outer_positions=None,
                                  antialiasing=True,
                                  colour_background=None))


class GUIMainWindow(QtGui.QMainWindow):

    def __init__(self):

        super(GUIMainWindow, self).__init__()

        self._image = None
        self.data_array = None
        self.set_loging(False)  # todo checkbox in settings
        self.settings = dialogs.SettingsDialog(self, DEFAULT_ARRAY[2])
        self.initUI()
        self.show()

    def set_loging(self, onoff):

        if onoff:
            self.logger = LogFile(log_filename="log/qt",
                                  override_log_files=True,
                                  log_colours=False,
                                  properties_different_colour=False)
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

        saveAction = QtGui.QAction('&Save stimulus', self)
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
        pil_generator = ICON[2]
        da_generator = ICON[1]
        self._image = pil_generator.plot(dot_array=da_generator.generate(n_dots=ICON[0]))

        self.setWindowIcon(QtGui.QIcon(self.pixmap()))
        self._image = None

        self.action_generate_btn()

    def make_new_array(self):

        try:
            generator = self.get_generator()
            self.data_array = generator.generate(n_dots=self.get_number(), logger=self.logger)
        except (RuntimeError, ValueError) as error:
            self.main_widget.text_error_feedback(error)
            raise error

        if self.settings.bicoloured.isChecked():
            data_array2 = generator.generate(n_dots=self.main_widget.number2.value,
                                             occupied_space=self.data_array,
                                             logger=self.logger)
            data_array2.attributes.change(colour=self.main_widget.dot_colour2.text)
            self.data_array.join(data_array2, realign=False)

        self._image = None

    def image(self):
        if self._image is not None:
            return self._image
        else:
            para = self.get_image_parameter()
            pil_gen = pil_image.PILImagePlotter(colour_target_area=para.colour_target_area,
                                                colour_field_area=para.colour_field_area,
                                                colour_field_area_outer=para.colour_field_area_outer,
                                                colour_center_of_mass=para.colour_center_of_mass,
                                                colour_center_of_outer_positions=para.colour_center_of_outer_positions,
                                                antialiasing=para.antialiasing,
                                                colour_background=para.colour_background)

            self._image = pil_gen.plot(dot_array=self.data_array)
            return self._image

    def get_number(self):
        return self.main_widget.number.value

    def get_generator(self):

        try:
            colour_dot = Colour(self.main_widget.dot_colour.text)
        except:
            colour_dot = DEFAULT_ARRAY[1].item_colour
            self.main_widget.dot_colour.text = colour_dot

        return DotArrayFactory(target_area_radius=self.main_widget.target_array_radius.value,
                               item_colour=colour_dot,
                               item_diameter_mean=self.main_widget.item_diameter_mean.value,
                               item_diameter_range=[self.main_widget.item_diameter_range.value1,
                                                      self.main_widget.item_diameter_range.value2],
                               item_diameter_std=self.main_widget.item_diameter_std.value,
                               minimum_gap=self.main_widget.minimum_gap.value)

    def get_image_parameter(self):
        # check colour input

        try:
            colour_area = Colour(self.settings.colour_area.text)
        except:
            colour_area = None
            self.settings.colour_area.text = "None"
        try:
            colour_convex_hull_positions = Colour(self.settings.colour_convex_hull_positions.text)
        except:
            colour_convex_hull_positions = None
            self.settings.colour_convex_hull_positions.text = "None"
        try:
            colour_convex_hull_dots = Colour(self.settings.colour_convex_hull_dots.text)
        except:
            colour_convex_hull_dots = None
            self.settings.colour_convex_hull_dots.text = "None"
        try:
            colour_background = Colour(self.settings.colour_background.text)
        except:
            colour_background = None
            self.settings.colour_background.text = "None"

        return pil_image.PILImagePlotter(colour_target_area=colour_area,
                                         colour_field_area=colour_convex_hull_positions,
                                         colour_field_area_outer=colour_convex_hull_dots,
                                         colour_center_of_mass=None,
                                         colour_center_of_outer_positions=None,
                                         antialiasing=self.settings.antialiasing.isChecked(),
                                         colour_background=colour_background)

    def pixmap(self):
        return QtGui.QPixmap.fromImage(ImageQt(self.image()))

    def show_current_image(self, remake_image=False):
        """"""
        if remake_image:
            self._image = None
        w = self.get_generator().target_array_radius * 2
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
        txt = self.data_array.get_features_text(extended_format=True, with_object_id=True)
        if self.settings.bicoloured.isChecked():
            prop = self.data_array.get_features_split_by_colours()
            for p in prop.split:
                txt += p.get_nice_text()
        if clear_field:
            self.main_widget.text_clear()

        self.main_widget.text_out(txt)

    def action_print_xy(self):
        """"""
        txt = self.data_array.get_csv(object_id_column=False, num_idx_column=False, colour_column=True)
        self.main_widget.text_out(txt)

    def action_print_para(self):
        d = {'number': self.get_number()}
        d['image_parameter'] = self.get_image_parameter().as_dict()
        d['generator'] = self.get_generator().as_dict()
        self.main_widget.text_out("# parameter\n" + yaml.dump(d, default_flow_style=False))

    def save_array(self):
        """"""

        filename, extension = QtGui.QFileDialog.getSaveFileNameAndFilter(
            self,
            'Save file',
            "",
            "Image PNG File (*.png);; Image BMP File (*.bmp);; JSON File (.json)"
        )

        if len(filename)>0:
            filename = os.path.abspath(filename)

            ext = extension.split("(")[1].replace(")", "")\
                                         .replace("*", "").strip()
            if not filename.endswith(ext):
                filename = filename + ext

            if ext == ".json":
                self.data_array.save(filename)

            elif ext == ".png" or ext == ".bmp":
                self.image().save(filename)


    def action_match(self):
        """"""
        prop = self.data_array.get_features_dict()
        response = dialogs.MatchPropertyDialog.get_response(self, prop)  #
        if response is not None:
            self.data_array.match(response, realign=True, center_array=True)
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
        self.data_array.attributes.change(colour=self.get_generator().item_attributes.colour)
        self.show_current_image(remake_image=True)

    def action_slider_released(self):
        """"""
        change = self.main_widget.number.value - self.data_array.feature_numerosity
        self.data_array = self.data_array.number_deviant(change)
        self.show_current_image(remake_image=True)
        self.write_properties()
        # todo slider does not work correctly for multi colour arrays

    def action_make_sequence(self):
        match_methods, match_range, extra_space = dialogs.SequenceDialog.get_response(self)

        d = {"match range": match_range, "extra_space": extra_space}
        d["match_methods"] = list(map(lambda x: x.as_dict(), match_methods))  # TODO <-- check ERROR
        self.main_widget.text_out("# Sequence\n" + \
                                           yaml.dump(d, default_flow_style=False))

        if match_methods is not None:
            # print("processing")
            sequence = generate_da_sequence(reference_dot_array=self.data_array,
                                   match_properties=match_methods,
                                   min_max_numerosity=match_range,
                                   extra_space=extra_space,
                                   logger=self.logger)
            SequenceDisplay(self, da_sequence=sequence,
                            start_numerosity=self.data_array.feature_numerosity,
                            image_parameter=self.get_image_parameter()).exec_()


