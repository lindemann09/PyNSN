"""
The named colour are the 140 HTML colour names:
   see https://www.w3schools.com/colors/colors_names.asp
"""
from functools import total_ordering

_NUMERALS = '0123456789abcdefABCDEF'
_HEXDEC = {v: int(v, 16) for v in (x + y for x in _NUMERALS for y in _NUMERALS)}


@total_ordering
class Colour(object):

    def __init__(self, colour=None):
        self.colour = colour

    def __repr__(self):
        return "Colour({})".format(self.colour)

    def __str__(self):
        return self._colour

    def __hash__(self):
        return hash(self._colour)

    def __lt__(self, other):
        return self._colour < other._colour

    def __eq__(self, other):
        return self._colour == other._colour

    def __ne__(self, other):
        return self._colour != other._colour

    @property
    def colour(self):
        return self._colour

    @colour.setter
    def colour(self, value):

        if value is None:
            self._colour = None
        elif isinstance(value, Colour):
            self._colour = value._colour
        else:
            try:
                Colour.hextriplet2rgb(value)  # check if vaild hextriplet
                if value[0] != "#":
                    value = "#" + value
                self._colour = value.upper()
            except:
                try:
                    self._colour = Colour.NAMED_COLOURS[value]
                except:
                    try:
                        self._colour = Colour.rgb2hextriplet(value)
                    except:
                        raise RuntimeError("Incorrect colour " + \
                                           "('{}')!\n Use RGB tuple, hex triplet or a colour name from Colour.NAMED_COLOURS.".format(
                                               value))

    @property
    def rgb(self):
        return Colour.hextriplet2rgb(self._colour)

    @staticmethod
    def hextriplet2rgb(hextriplet):
        ht = hextriplet.lstrip("#")
        return _HEXDEC[ht[0:2]], _HEXDEC[ht[2:4]], _HEXDEC[ht[4:6]]

    @staticmethod
    def rgb2hextriplet(rgb, uppercase=True):
        if len(rgb) != 3:
            raise TypeError("rgb must be a list of three values.")
        if uppercase:
            lettercase = 'X'
        else:
            lettercase = 'x'
        return "#" + format(rgb[0] << 16 | rgb[1] << 8 | rgb[2], '06' + lettercase)

    NAMED_COLOURS = {

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

    def __init__(self,
                 target_area=None,
                 field_area=None,
                 field_area_outer=None,
                 center_of_mass=None,
                 center_of_outer_positions=None,
                 background=None,
                 default_dot_colour=Colour("lightgreen"),  # used if no color
                 # specified in dot array
                 ):

        self.target_area = Colour(target_area)
        self.field_area = Colour(field_area)
        self.field_area_outer = Colour(field_area_outer)
        self.center_of_mass = Colour(center_of_mass)
        self.center_of_outer_positions = Colour(center_of_outer_positions)
        self.background = Colour(background)
        self.default_dot_colour = Colour(default_dot_colour)

    def as_dict(self):
        return {"colour_total_area": self.target_area.colour,
                "colour_field_area": self.field_area.colour,
                "colour_field_area_outer": self.field_area_outer.colour,
                "colour_center_of_mass": self.center_of_mass.colour,
                "colour_center_of_outer_positions": self.center_of_outer_positions.colour,
                "colour_background": self.background.colour,
                "default_dot_colour": self.default_dot_colour.colour}

