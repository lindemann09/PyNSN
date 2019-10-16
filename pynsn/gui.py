#!/usr/bin/python

"""
"""
try:
    from PyQt4 import QtGui
except:
    raise ModuleNotFoundError("Running the PyNSN GUI requires the installation 'PyQt4' (version 4.10 or larger)")

import sys
from .qt.gui_main_window import GUIMainWindow

def start():
    app = QtGui.QApplication(sys.argv)
    ex = GUIMainWindow()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    start()
