
import json
import os

from PIL.ImageQt import ImageQt
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtWidgets import QAction, QMainWindow, qApp, QFileDialog

from . import dialogs
from .main_widget import MainWidget
from .sequence_display import SequenceDisplay
from .. import __version__
from .. import random_array
from .. import arrays
from .._lib.exceptions import NoSolutionError
from ..random_array import distributions as distr
from .. import visual_properties as props

from ..image import _colour
from ..image import pil_image
from .._sequence import dot_array_sequence

DEFAULT_ARRAY = {"num": 40,
        "ref": arrays.BaseArray(target_area_radius=200,
                                min_dist_between=2),
        "sdr": random_array.SizeDistribution(
                    dot_diameter=distr.Beta(mu=15, sigma=8, min_max=(5, 40))),
        "col": _colour.ImageColours(target_area="#303030",
                                    field_area_positions=None,
                                    field_area=None,
                                    center_of_field_area=None,
                                    center_of_mass=None,
                                    default_object_colour="green",
                                    background="gray")}

ICON = {"num": 11,
        "ref": arrays.BaseArray(target_area_radius=200),
        "sdr": random_array.SizeDistribution(
                dot_diameter=distr.Beta(mu=35, sigma=20, min_max=(5, 80))),
        "col": _colour.ImageColours(target_area="#3e3e3e",
                                    field_area_positions=None,
                                    field_area="expyriment_orange",
                                    center_of_field_area=None,
                                    center_of_mass=None,
                                    default_object_colour="lime",
                                    background=None)}


class GUIMainWindow(QMainWindow):

    def __init__(self):

        super(GUIMainWindow, self).__init__()

        self._image = None
        self.dot_array = None
        self.settings = dialogs.SettingsDialog(self, image_colours=DEFAULT_ARRAY["col"])
        self.initUI()
        self.show()


    def initUI(self):

        # menus
        exitAction = QAction(QIcon('exit.png'), '&Exit',
                                       self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(qApp.quit)

        settingsAction = QAction('&Settings', self)
        settingsAction.triggered.connect(self.action_settings)

        saveAction = QAction('&Save current image', self)
        saveAction.setShortcut('Ctrl+S')
        saveAction.triggered.connect(self.save_array)


        printxyAction = QAction('&Print array', self)
        printxyAction.triggered.connect(self.action_print_xy)
        printparaAction = QAction('&Print parameter', self)
        printparaAction.triggered.connect(self.action_print_para)

        adaptAction = QAction('&Adapt property', self)
        adaptAction.triggered.connect(self.action_adapt)

        sequenceAction = QAction('&Make _sequence', self)
        sequenceAction.triggered.connect(self.action_make_sequence)

        aboutAction = QAction('&About', self)

        # self.statusBar()
        menubar = self.menuBar()

        arrayMenu = menubar.addMenu('&Array')
        arrayMenu.addAction(settingsAction)
        arrayMenu.addAction(saveAction)
        arrayMenu.addSeparator()
        arrayMenu.addAction(exitAction)

        toolMenu = menubar.addMenu('&Tools')
        toolMenu.addAction(sequenceAction)
        toolMenu.addAction(adaptAction)
        toolMenu.addSeparator()
        toolMenu.addAction(printparaAction)
        toolMenu.addAction(printxyAction)

        aboutMenu = menubar.addMenu('&About')
        aboutMenu.addAction(aboutAction)



        # main widget
        self.main_widget = MainWidget(self, self.settings, DEFAULT_ARRAY["num"],
                                      DEFAULT_ARRAY["ref"], DEFAULT_ARRAY["sdr"])
        self.setCentralWidget(self.main_widget)
        self.main_widget.btn_generate.clicked.connect(self.action_generate_btn)
        self.main_widget.dot_colour.edit.editingFinished.connect(self.action_dot_colour_change)

        self.move(300, -300)
        self.setWindowTitle('PyNSN GUI {}'.format(__version__))

        # ICON
        oa = random_array.create(reference_array=ICON["ref"],
                                 size_distribution=ICON["sdr"],
                                 n_objects=ICON["num"])
        self._image = pil_image.create( object_array=oa,
                        colours=ICON["col"], antialiasing=True)

        self.setWindowIcon(QIcon(self.pixmap()))
        self._image = None
        self.action_generate_btn()

    def make_new_array(self):
        ref_ar, sdr = self.get_specs()
        try:
            self.dot_array = random_array.create(n_objects=self.get_number(),
                                                 reference_array=ref_ar,
                                                 size_distribution=sdr)
        except NoSolutionError as error:
            self.main_widget.text_error_feedback(error)
            raise error

        if self.settings.bicoloured.isChecked():
            data_array2 = random_array.create(n_objects=self.main_widget.number2.value,
                                              reference_array=ref_ar,
                                              size_distribution=sdr,
                                              occupied_space=self.dot_array)
            data_array2.set_attributes(self.main_widget.dot_colour2.text)
            self.dot_array.join(data_array2)

        self.dot_array.round(decimals=self.settings.rounding_decimals.value)
        self._image = None

    def image(self):
        if self._image is not None:
            return self._image
        else:
            para = self.get_image_colours()
            image_colours = _colour.ImageColours(
                target_area=para.target_area,
                field_area_positions=para.field_area_positions,
                field_area=para.field_area,
                center_of_field_area=para.center_of_field_area,
                center_of_mass=para.center_of_mass,
                default_object_colour=para.default_object_colour,
                background=para.background)

            self._image = pil_image.create(object_array=self.dot_array,
                                           colours=image_colours,
                                           antialiasing=self.settings.antialiasing.isChecked())
                                            # todo maybe: gabor_filter=ImageFilter.GaussianBlur

            return self._image

    def get_number(self):
        return self.main_widget.number.value

    def get_specs(self):
        size_dist = random_array.SizeDistribution(dot_diameter=
                        distr.Beta(mu=self.main_widget.item_diameter_mean.value,
                       sigma=self.main_widget.item_diameter_std.value,
                       min_max=[self.main_widget.item_diameter_range.value1,
                                self.main_widget.item_diameter_range.value2])
                                                  )
        ref_array = arrays.BaseArray(
            target_area_radius=self.main_widget.target_area_radius.value,
            min_dist_between=self.main_widget.min_dist_between.value)
        return ref_array, size_dist

    def get_image_colours(self):
        # check colour input
        try:
            colour_area = _colour.Colour(self.settings.colour_area.text)
        except TypeError:
            colour_area = None
            self.settings.colour_area.text = "None"
        try:
            colour_convex_hull_positions = _colour.Colour(
                self.settings.colour_convex_hull_positions.text)
        except TypeError:
            colour_convex_hull_positions = None
            self.settings.colour_convex_hull_positions.text = "None"
        try:
            colour_convex_hull_dots = _colour.Colour(
                self.settings.colour_convex_hull_dots.text)
        except TypeError:
            colour_convex_hull_dots = None
            self.settings.colour_convex_hull_dots.text = "None"
        try:
            colour_background = _colour.Colour(
                self.settings.colour_background.text)
        except TypeError:
            colour_background = None
            self.settings.colour_background.text = "None"

        return _colour.ImageColours(target_area=colour_area,
                                    field_area_positions=colour_convex_hull_positions,
                                    field_area=colour_convex_hull_dots,
                                    center_of_field_area=None,
                                    center_of_mass=None,
                                    default_object_colour=self.settings.default_object_colour,
                                    background=colour_background)

    def pixmap(self):
        return QPixmap.fromImage(ImageQt(self.image()))

    def show_current_image(self, remake_image=False):
        """"""
        if remake_image:
            self._image = None
        ref_ar, _ = self.get_specs()
        w = ref_ar.target_area_radius * 2
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
        txt = self.dot_array._properties.as_text(extended_format=True,
                                               with_hash=True)
        if self.settings.bicoloured.isChecked():
            for da in self.dot_array.get_split_arrays():
                txt += "Attribute {}\n".format(da._attributes[0])
                txt += da._properties.as_text(extended_format=True,
                                            with_hash=False)
        if clear_field:
            self.main_widget.text_clear()

        self.main_widget.text_out(txt)

    def action_print_xy(self):
        """"""
        txt = self.dot_array.csv(hash_column=False, num_idx_column=False, colour_column=True)
        self.main_widget.text_out(txt)

    def action_print_para(self):
        d = {'number': self.get_number()}
        d['image_parameter'] = self.get_image_colours().as_dict()
        ref_array, size_distr = self.get_specs()
        d['reference_array'] = ref_array.as_dict()
        d['size_distribution'] = size_distr.as_dict()
        self.main_widget.text_out("# parameter\n" + json.dumps(d, indent=2))

    def save_array(self):
        """"""

        filename, extension = QFileDialog.getSaveFileName(
            self,
            caption='Save file',
            directory=".",
            filter="Image PNG File (*.png);; Image BMP File (*.bmp);; JSON File (.json)")

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

    def action_adapt(self):
        """"""
        prop = self.dot_array._properties.as_dict()
        prop, value = dialogs.AdaptPropertyDialog.get_response(self, prop)
        if prop is not None:
            self.dot_array = props.fit.visual_property(self.dot_array,
                                                prop, value=value)
            if prop in [x for x in props.flags if x.is_size_property()]:
                self.dot_array.center_array()
                self.dot_array.realign()
            self.show_current_image(remake_image=True)
            self.write_properties()
            self.main_widget.updateUI()

    def action_settings(self):
        """"""
        result = self.settings.exec_()
        self.main_widget.updateUI()
        self.action_generate_btn()
        #self.show_current_image(remake_image=True)

    def action_dot_colour_change(self):
        """
        """
        try:
            colour_dot = _colour.Colour(self.main_widget.dot_colour.text)
            self.settings.default_object_colour = colour_dot
        except TypeError:
            colour_dot = self.settings.default_object_colour
            self.main_widget.dot_colour.text = colour_dot.colour

        self.dot_array.set_attributes(colour_dot.colour)
        self.show_current_image(remake_image=True)

    def action_make_sequence(self):
        adapt_methods, adapt_range, extra_space = \
                                    dialogs.SequenceDialog.get_response(self)
        adapt_methods = adapt_methods[0] #FIXME just adapt the first

        d = {"adapt range": adapt_range,
             "extra_space": extra_space}
        d["adapt_methods"] = adapt_methods
        self.main_widget.text_out("# Sequence\n" + \
                                           json.dumps(d))
        ref_array, sdr = self.get_specs()
        ref_array.min_dist_area_boarder = extra_space/2
        ref_array.target_area_radius += ref_array.min_dist_area_boarder

        if adapt_methods is not None:
            sequence = dot_array_sequence.create(
                              specs=specs,
                              adapt_property=adapt_methods,
                              adapt_value=self.dot_array._properties.get(adapt_methods),
                              min_max_numerosity=adapt_range)
            SequenceDisplay(self, da_sequence=sequence,
                            start_numerosity=self.dot_array._properties.numerosity,
                            image_colours=self.get_image_colours(),
                            antialiasing=self.settings.antialiasing.isChecked()).exec_()


