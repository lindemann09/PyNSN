"""Creating Non-Symbolic Number Displays"""
# pylint: disable=C0413

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
__version__ = '0.9.dev1'

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

from ._shapes import (Point2D, Dot, Rectangle, Picture, Ellipse, PolygonShape,
                      Colour)
from ._stimulus import VisProp, NSNStimulus, StimulusFactory
