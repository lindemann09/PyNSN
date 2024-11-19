"""Creating Non-Symbolic Number Displays"""
# pylint: disable=C0413

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
__version__ = '1.0.6'

from sys import version_info as _python_version_info
from ._misc import is_interactive_mode as _is_interactive_mode

if not (_python_version_info[0] >= 3 and _python_version_info[1] >= 10):
    raise RuntimeError(
        f"PyNSN {__version__} is not compatible with Python " +
        f"{_python_version_info[0]}.{_python_version_info[1]}. " +
        "Please use Python 3.10 or later."
    )
if _is_interactive_mode():
    print(f"PyNSN {__version__}")

from . import defaults
from . import exceptions
from ._shapes import (Point2D, Dot, Rectangle, Picture, Ellipse, PolygonShape,
                      Colour)
from ._stimulus import (NSNStimulus, NSNStimulusPair, VisProp,
                        SPACE_PROPERTIES, SIZE_PROPERTIES)
from ._factory import StimulusFactory
from . import rnd
from . import fit
