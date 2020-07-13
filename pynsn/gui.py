#!/usr/bin/python

"""
"""
try:
    from PyQt5.QtWidgets import QApplication as _QApplication
except:
    raise ModuleNotFoundError("Running the PyNSN GUI requires the "
                              "installation 'PyQt5'")

import sys as _sys

from ._qt.gui_main_window import GUIMainWindow as _GUIMainWindow


def start():
    app = _QApplication(_sys.argv)
    ex = _GUIMainWindow()
    ex.show()
    _sys.exit(app.exec_())
