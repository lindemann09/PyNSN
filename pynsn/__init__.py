"""Creating Non-Symbolic Number Displays"""

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
__version__ = "1.1.1-dev2"

__all__ = (
    "defaults",
    "exceptions",
    "Point2D",
    "Dot",
    "Rectangle",
    "Picture",
    "Ellipse",
    "PolygonShape",
    "Colour",
    "NSNStimulus",
    "NSNStimulusPair",
    "VP",
    "rnd",
    "StimulusFactory",
    "fit",
)

from sys import version_info as _python_version_info

from . import defaults, exceptions, fit, rnd
from ._misc import is_interactive_mode as _is_interactive_mode
from ._shapes import Colour, Dot, Ellipse, Picture, Point2D, PolygonShape, Rectangle
from ._stimulus import VP, NSNStimulus, NSNStimulusPair
from .rnd._factory import StimulusFactory


if not (_python_version_info[0] >= 3 and _python_version_info[1] >= 10):
    raise RuntimeError(
        f"PyNSN {__version__} is not compatible with Python "
        + f"{_python_version_info[0]}.{_python_version_info[1]}. "
        + "Please use Python 3.10 or later."
    )
if _is_interactive_mode():
    print(f"PyNSN {__version__}")
