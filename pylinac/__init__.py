import sys

from pylinac.version import __version__

# check python version
if sys.version_info[0] < 3 or sys.version_info[1] < 7:
    raise ValueError(
        "Pylinac is only supported on Python 3.7+. Please update your environment."
    )

# import shortcuts
# core first
from .core import decorators, geometry, image, io, mask, profile, roi, utilities

# alphabetized modules
from .acr import ACRCT, ACRMRILarge
from .ct import CatPhan504, CatPhan600, CatPhan503, CatPhan604
from .quart import QuartDVT
from .core.utilities import clear_data_files, assign2machine
from .field_analysis import (
    FieldAnalysis,
    DeviceFieldAnalysis,
    Protocol,
    Device,
    Edge,
    Interpolation,
    Normalization,
    Centering,
)
from .planar_imaging import (
    LeedsTOR,
    StandardImagingQC3,
    LasVegas,
    DoselabMC2kV,
    DoselabMC2MV,
    StandardImagingQCkV,
    PTWEPIDQC,
    SNCMV,
    SNCkV,
    StandardImagingFC2,
    IMTLRad,
    SNCFSQA,
    LeedsTORBlue,
)
from .log_analyzer import load_log, Dynalog, TrajectoryLog, MachineLogs
from .picketfence import PicketFence  # must be after log analyzer
from .starshot import Starshot
from .vmat import DRMLC, DRGS
from .winston_lutz import WinstonLutz
from .calibration import tg51, trs398
