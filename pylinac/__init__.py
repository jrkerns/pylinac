
import sys

__version__ = '2.5.0.0'
__version_info__ = (2, 5, 0, 0)

# check python version
if sys.version_info[0] < 3 or sys.version_info[1] < 6:
    raise ValueError("Pylinac is only supported on Python 3.6+. Please update your environment.")

# import shortcuts
from pylinac.ct import CatPhan500, CatPhan504, CatPhan600, CatPhan503, CatPhan604
from pylinac.core import decorators, geometry, image, io, mask, profile, roi, utilities
from pylinac.core.utilities import clear_data_files, assign2machine
from pylinac.flatsym import FlatSym
from pylinac.planar_imaging import LeedsTOR, StandardImagingQC3, LasVegas, DoselabMC2kV, DoselabMC2MV
from pylinac.log_analyzer import load_log, Dynalog, TrajectoryLog, MachineLogs
from pylinac.picketfence import PicketFence  # must be after log analyzer
from pylinac.starshot import Starshot
from pylinac.vmat import DRMLC, DRGS
from pylinac.winston_lutz import WinstonLutz
from pylinac.calibration import tg51, trs398

