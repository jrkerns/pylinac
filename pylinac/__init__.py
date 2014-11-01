

from __future__ import absolute_import
import sys

__version__ = '0.1.3'
__version_info__ = (0,1,3)

# Determine if using python 2 or 3 (mostly for incompatible name clashes like Tkinter/tkinter)
running_py3 = sys.version_info[0] == 3

# Determine if user has PySide. If they do, Dialogs will use that; if not, will use tkinter
try:
    import PyQt4
    has_pyside = True
except ImportError:
    has_pyside = False

# import major packages into main namespace for convenience
from pylinac.starshot.starshot import Starshot
from pylinac.vmatqa.vmat import VMAT



