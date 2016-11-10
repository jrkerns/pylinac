
import sys

__version__ = '1.7.2'
__version_info__ = (1, 7, 2)

# check for python 2
if sys.version_info[0] < 3:
    raise ValueError("Pylinac is only supported on Python 3.x. It seems you are using Python 2; please use a different interpreter.")

# import shortcuts
from pylinac.cbct import CBCT
from pylinac.core import decorators, geometry, image, io, mask, profile, roi, utilities
from pylinac.core.utilities import clear_data_files
from pylinac.flatsym import BeamImage
from pylinac.planar_imaging import LeedsTOR, StandardImagingQC3
from pylinac.log_analyzer import load_log, Dynalog, TrajectoryLog, MachineLogs
from pylinac.picketfence import PicketFence
from pylinac.starshot import Starshot
from pylinac.vmat import VMAT
from pylinac.winston_lutz import WinstonLutz
from pylinac import tg51

from pylinac.watcher import watch, process
