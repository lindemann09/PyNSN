"""PyNSN package

Creating Non-Symbolic Number Displays
"""

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"
__version__ = "1.0-dev"

from sys import version_info as _python_version_info

# pylint: disable=C0209
if not (_python_version_info[0] >= 3 and _python_version_info[1] >= 8):
    raise RuntimeError(
        "PyNSN {0} ".format(__version__)
        + "is not compatible with Python {0}.{1}. ".format(
            _python_version_info[0], _python_version_info[1]
        )
        + "Please use Python 3.8 or later."
    )

from ._shapes import (Point2D, Dot, Rectangle, Picture,
                      Ellipse, PolygonShape, ShapeType, CircularShapeType, Colour)
from ._stimulus import ShapeArray, ArrayProperties, VisProp, NSNStimulus


def _print_version_info():
    from ._misc import is_interactive_mode

    if is_interactive_mode():
        print(f"PyNSN {__version__}")


_print_version_info()
