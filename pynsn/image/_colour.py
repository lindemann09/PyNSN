"""
The named colour are the 140 HTML colour names:
   see https://www.w3schools.com/colors/colors_names.asp
"""
from collections import OrderedDict
from functools import total_ordering
from typing import Sequence, Optional, Union, Tuple

from .._lib import misc

_NUMERALS = '0123456789abcdefABCDEF'
_HEXDEC = {v: int(v, 16)
           for v in (x + y for x in _NUMERALS for y in _NUMERALS)}


RGBType = Tuple[int, int, int]
ColourType = Union[str, RGBType, Sequence[float], None]


@ total_ordering
class Colour(object):
    """Colour Class

    Args:
        colour:
            Hextriplet, RGB value, colour name and ``Colour``.

            If colour is unknown string or ``None``, the specified default colour will be used.
            Thus, ``Colour(variable, default="red")`` will result in 'red', if `variable`
            is unknown or ``None``. If default is not defined the colour will be a
            'None Colour'.

        default:
            default colour


    Notes:
        The following colour names are known:

        .. code-block:: py

            'aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'black',
            'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse',
            'chocolate', 'coral', 'cornflowerblue', 'cornsilk', 'crimson', 'cyan', 'darkblue',
            'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkkhaki', 'darkmagenta',
            'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen',
            'darkslateblue', 'darkslategray', 'darkturquoise', 'darkviolet', 'deeppink',
            'deepskyblue', 'dimgray', 'dodgerblue', 'firebrick', 'floralwhite', 'forestgreen',
            'fuchsia', 'gainsboro', 'ghostwhite', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow',
            'honeydew', 'hotpink', 'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush',
            'lawngreen', 'lemonchiffon', 'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow',
            'lightgreen', 'lightgray', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue',
            'lightslategray', 'lightsteelblue', 'lightyellow', 'lime', 'limegreen', 'linen', 'magenta',
            'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen',
            'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue',
            'mintcream', 'mistyrose', 'moccasin', 'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab',
            'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred',
            'papayawhip', 'peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple', 'red', 'rosybrown',
            'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', 'silver',
            'skyblue', 'slateblue', 'slategray', 'snow', 'springgreen', 'steelblue', 'tan', 'teal',
            'thistle', 'tomato', 'turquoise', 'violet', 'wheat', 'white', 'whitesmoke', 'yellow',
            'yellowgreen', 'expyriment_orange', 'expyriment_purple'
    """

    def __init__(self,
                 colour: ColourType,
                 default: ColourType = None) -> None:
        if colour is None:
            colour = default

        self._colour = None
        try:
            self.set(colour)
        except TypeError:  # colour is not a valid colour
            self.set(default)

    def __repr__(self) -> str:
        return "Colour({})".format(self.colour)

    def __hash__(self) -> int:
        return hash(self._colour)

    def __lt__(self, other):
        return self._colour < other._colour

    def __eq__(self, other):
        return self._colour == other._colour

    def __ne__(self, other):
        return self._colour != other._colour

    @property
    def colour(self) -> Union[str, None]:
        """Hextriplet code of the colour or None"""
        return self._colour

    def set(self, value: ColourType) -> None:
        """Set the colour values

        Args:
            value:

        Returns:

        """
        error_txt = "Incorrect colour ('{}')!\n Use RGB tuple, " +\
                    "hex triplet or a colour name from Colour.NAMED_COLOURS."

        if value is None:
            self._colour = None
        elif isinstance(value, Colour):
            self._colour = value.colour

        elif isinstance(value, str):
            try:
                # check if valid hextriplet
                Colour.convert_hextriplet_to_rgb(value)
                if value[0] != "#":
                    value = "#" + value
                self._colour = value.upper()
            except ValueError:
                try:
                    self._colour = Colour.NAMED_COLOURS[value]
                except KeyError as err:
                    raise TypeError(error_txt.format(value)) from err
        else:
            try:
                self._colour = Colour.convert_rgb_hextriplet(value)
            except TypeError as err:
                raise TypeError(error_txt.format(value)) from err

    @property
    def rgb(self) -> Union[RGBType, None]:
        """RGB code of the colour"""
        if self._colour is not None:
            return Colour.convert_hextriplet_to_rgb(self._colour)

    def get_rgb_alpha(self, alpha: float) -> Union[Tuple[int, int, int, int], None]:
        """RBG with alpha values

        Args:
            alpha:
             Alpha can be float <= 1.0 or an integer between [0, 255]

        Returns:
            tuple of four values, RGB and alpha value
        """

        if self.rgb is None:
            return None

        if isinstance(alpha, float):
            if alpha < 0. or alpha > 1.:
                raise TypeError(
                    "If alpha is a float, it has to be between 0 and 1.")
            alpha = round(alpha * 255)
        elif isinstance(alpha, int):
            if alpha < 0 or alpha > 255:
                raise TypeError(
                    "If alpha is an int, it has to be between 0 and 255.")
        else:
            raise TypeError("If alpha has to be a float or int and not {}.".format(
                type(alpha)))

        return self.rgb + (alpha,)

    @ staticmethod
    def convert_hextriplet_to_rgb(hextriplet: str) -> RGBType:
        """Convert a hextriplet string to RGB values"""
        ht = hextriplet.lstrip("#")
        try:
            return _HEXDEC[ht[0:2]], _HEXDEC[ht[2:4]], _HEXDEC[ht[4:6]]
        except KeyError as err:
            raise ValueError(
                f"Can't convert {hextriplet} to rgb value") from err

    @ staticmethod
    def convert_rgb_hextriplet(rgb: Sequence, uppercase: bool = True):
        """Convert RBG values to a hextriplet string"""
        if len(rgb) != 3:
            raise TypeError("rgb must be a list of three values.")
        if uppercase:
            lettercase = 'X'
        else:
            lettercase = 'x'

        return "#" + format(rgb[0] << 16 | rgb[1] << 8 | rgb[2], '06' + lettercase)

    NAMED_COLOURS = {  # Dict with known colour names and hextriplets
        'aliceblue': '#F0F8FF',
        'antiquewhite': '#FAEBD7',
        'aqua': '#00FFFF',
        'aquamarine': '#7FFFD4',
        'azure': '#F0FFFF',
        'beige': '#F5F5DC',
        'bisque': '#FFE4C4',
        'black': '#000000',
        'blanchedalmond': '#FFEBCD',
        'blue': '#0000FF',
        'blueviolet': '#8A2BE2',
        'brown': '#A52A2A',
        'burlywood': '#DEB887',
        'cadetblue': '#5F9EA0',
        'chartreuse': '#7FFF00',
        'chocolate': '#D2691E',
        'coral': '#FF7F50',
        'cornflowerblue': '#6495ED',
        'cornsilk': '#FFF8DC',
        'crimson': '#DC143C',
        'cyan': '#00FFFF',
        'darkblue': '#00008B',
        'darkcyan': '#008B8B',
        'darkgoldenrod': '#B8860B',
        'darkgray': '#A9A9A9',
        'darkgreen': '#006400',
        'darkkhaki': '#BDB76B',
        'darkmagenta': '#8B008B',
        'darkolivegreen': '#556B2F',
        'darkorange': '#FF8C00',
        'darkorchid': '#9932CC',
        'darkred': '#8B0000',
        'darksalmon': '#E9967A',
        'darkseagreen': '#8FBC8F',
        'darkslateblue': '#483D8B',
        'darkslategray': '#2F4F4F',
        'darkturquoise': '#00CED1',
        'darkviolet': '#9400D3',
        'deeppink': '#FF1493',
        'deepskyblue': '#00BFFF',
        'dimgray': '#696969',
        'dodgerblue': '#1E90FF',
        'firebrick': '#B22222',
        'floralwhite': '#FFFAF0',
        'forestgreen': '#228B22',
        'fuchsia': '#FF00FF',
        'gainsboro': '#DCDCDC',
        'ghostwhite': '#F8F8FF',
        'gold': '#FFD700',
        'goldenrod': '#DAA520',
        'gray': '#808080',
        'green': '#008000',
        'greenyellow': '#ADFF2F',
        'honeydew': '#F0FFF0',
        'hotpink': '#FF69B4',
        'indianred': '#CD5C5C',
        'indigo': '#4B0082',
        'ivory': '#FFFFF0',
        'khaki': '#F0E68C',
        'lavender': '#E6E6FA',
        'lavenderblush': '#FFF0F5',
        'lawngreen': '#7CFC00',
        'lemonchiffon': '#FFFACD',
        'lightblue': '#ADD8E6',
        'lightcoral': '#F08080',
        'lightcyan': '#E0FFFF',
        'lightgoldenrodyellow': '#FAFAD2',
        'lightgreen': '#90EE90',
        'lightgray': '#D3D3D3',
        'lightpink': '#FFB6C1',
        'lightsalmon': '#FFA07A',
        'lightseagreen': '#20B2AA',
        'lightskyblue': '#87CEFA',
        'lightslategray': '#778899',
        'lightsteelblue': '#B0C4DE',
        'lightyellow': '#FFFFE0',
        'lime': '#00FF00',
        'limegreen': '#32CD32',
        'linen': '#FAF0E6',
        'magenta': '#FF00FF',
        'maroon': '#800000',
        'mediumaquamarine': '#66CDAA',
        'mediumblue': '#0000CD',
        'mediumorchid': '#BA55D3',
        'mediumpurple': '#9370DB',
        'mediumseagreen': '#3CB371',
        'mediumslateblue': '#7B68EE',
        'mediumspringgreen': '#00FA9A',
        'mediumturquoise': '#48D1CC',
        'mediumvioletred': '#C71585',
        'midnightblue': '#191970',
        'mintcream': '#F5FFFA',
        'mistyrose': '#FFE4E1',
        'moccasin': '#FFE4B5',
        'navajowhite': '#FFDEAD',
        'navy': '#000080',
        'oldlace': '#FDF5E6',
        'olive': '#808000',
        'olivedrab': '#6B8E23',
        'orange': '#FFA500',
        'orangered': '#FF4500',
        'orchid': '#DA70D6',
        'palegoldenrod': '#EEE8AA',
        'palegreen': '#98FB98',
        'paleturquoise': '#AFEEEE',
        'palevioletred': '#DB7093',
        'papayawhip': '#FFEFD5',
        'peachpuff': '#FFDAB9',
        'peru': '#CD853F',
        'pink': '#FFC0CB',
        'plum': '#DDA0DD',
        'powderblue': '#B0E0E6',
        'purple': '#800080',
        'red': '#FF0000',
        'rosybrown': '#BC8F8F',
        'royalblue': '#4169E1',
        'saddlebrown': '#8B4513',
        'salmon': '#FA8072',
        'sandybrown': '#FAA460',
        'seagreen': '#2E8B57',
        'seashell': '#FFF5EE',
        'sienna': '#A0522D',
        'silver': '#C0C0C0',
        'skyblue': '#87CEEB',
        'slateblue': '#6A5ACD',
        'slategray': '#708090',
        'snow': '#FFFAFA',
        'springgreen': '#00FF7F',
        'steelblue': '#4682B4',
        'tan': '#D2B48C',
        'teal': '#008080',
        'thistle': '#D8BFD8',
        'tomato': '#FF6347',
        'turquoise': '#40E0D0',
        'violet': '#EE82EE',
        'wheat': '#F5DEB3',
        'white': '#FFFFFF',
        'whitesmoke': '#F5F5F5',
        'yellow': '#FFFF00',
        'yellowgreen': '#9ACD32',
        'expyriment_orange': '#FF9632',
        'expyriment_purple': '#A046FA'
    }


class ImageColours(object):
    COL_TARGET_AREA = "#FFF0D9"
    COL_DEFAULT_OBJECT = "darkgreen"
    OPACITY_OBJECT = 1
    OPACITY_GUIDES = 0.5

    def __init__(self,
                 target_area: ColourType = None,
                 field_area_positions: ColourType = None,
                 field_area: ColourType = None,
                 center_of_field_area: ColourType = None,
                 center_of_mass: ColourType = None,
                 background: ColourType = None,
                 default_object_colour: ColourType = None,
                 opacity_object: Optional[float] = None,
                 opacity_guides: Optional[float] = None
                 ):

        self.target_area = Colour(target_area,
                                  default=ImageColours.COL_TARGET_AREA)
        self.field_area_positions = Colour(field_area_positions, default=None)
        self.field_area = Colour(field_area, default=None)
        self.center_of_field_area = Colour(center_of_field_area, default=None)
        self.center_of_mass = Colour(center_of_mass, default=None)
        self.background = Colour(background, default=None)
        self.default_object_colour = Colour(default_object_colour,
                                            default=ImageColours.COL_DEFAULT_OBJECT)
        if opacity_guides is None:
            opacity_guides = ImageColours.OPACITY_GUIDES
        if opacity_guides < 0 or opacity_guides > 1:
            raise ValueError(
                f"opacity_guides ({opacity_guides}) has to be between 0 and 1")
        if opacity_object is None:
            opacity_object = ImageColours.OPACITY_OBJECT
        if opacity_object < 0 or opacity_object > 1:
            raise ValueError(
                f"opacity_object ({opacity_object}) has to be between 0 and 1")
        self.opacity_object = opacity_object
        self.opacity_guides = opacity_guides

    def to_dict(self) -> dict:
        return OrderedDict(
            {"total_area": self.target_area.colour,
             "field_area_positions": self.field_area_positions.colour,
             "field_area": self.field_area.colour,
             "center_of_field_area": self.center_of_field_area.colour,
             "center_of_mass": self.center_of_mass.colour,
             "background": self.background.colour,
             "default_object": self.default_object_colour.colour,
             "object_opacity": self.opacity_object,
             "info_opacity": self.opacity_guides
             })

    def __str__(self) -> str:
        return misc.dict_to_text(self.to_dict(), col_a=24, col_b=10)
