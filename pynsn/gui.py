"""
"""
try:
    from PyQt5.QtWidgets import QApplication as _QApplication
except:
    print("""ERROR: Running the PyNSN GUI requires the installation `PyQt5`.
Install `PyQt5` or  reinstall PyNSN via: `pip install pynsn[gui] -U`""")
    exit()

import sys as _sys
from ._gui.gui_main_window import GUIMainWindow as _GUIMainWindow

def start():
    app = _QApplication(_sys.argv)
    ex = _GUIMainWindow()
    ex.show()
    _sys.exit(app.exec_())

