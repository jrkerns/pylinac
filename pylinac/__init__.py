

from __future__ import absolute_import
import sys

__version__ = '0.1.0'
__version_info__ = (0,1,0)

# Determine if using python 2 or 3 (mostly for incompatible name clashes like Tkinter/tkinter)
running_py3 = sys.version_info[0] == 3

# Determine if user has PyQt4. If they do, Dialogs will use that; if not, will use tkinter
try:
    import PyQt4
    has_pyqt = True
except ImportError:
    has_pyqt = False

# import major packages into main namespace for convenience
from pylinac.starshot.starshot import Starshot
from pylinac.vmatqa.vmat import VMAT



