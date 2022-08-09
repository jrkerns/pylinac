
import sys

__version__ = '3.1.0'
__version_info__ = (3, 1, 0)

# check python version
if sys.version_info[0] < 3 or sys.version_info[1] < 6:
    raise ValueError("Pylinac is only supported on Python 3.6+. Please update your environment.")

# import shortcuts
from .ct import CatPhan504, CatPhan600, CatPhan503, CatPhan604
from .core import decorators, geometry, image, io, mask, profile, roi, utilities
from .core.utilities import clear_data_files, assign2machine
from .flatsym import FlatSym
from .planar_imaging import LeedsTOR, StandardImagingQC3, LasVegas, DoselabMC2kV, DoselabMC2MV
from .log_analyzer import load_log, Dynalog, TrajectoryLog, MachineLogs
from .picketfence import PicketFence  # must be after log analyzer
from .starshot import Starshot
from .vmat import DRMLC, DRGS
from .winston_lutz import WinstonLutz
from . import calibration
