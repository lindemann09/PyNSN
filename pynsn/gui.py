#!/usr/bin/python

"""
"""
try:
    from PyQt4 import QtGui as _QtGui
except:
    raise ModuleNotFoundError("Running the PyNSN GUI requires the installation 'PyQt4' (version 4.10 or larger)")

import sys as _sys
from ._qt.gui_main_window import GUIMainWindow as _GUIMainWindow

def start():
    app = _QtGui.QApplication(_sys.argv)
    ex = _GUIMainWindow()
    ex.show()
    _sys.exit(app.exec_())
