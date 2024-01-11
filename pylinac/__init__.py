import sys

from pylinac.version import __version__

# check python version
if sys.version_info[0] < 3 or sys.version_info[1] < 7:
    raise ValueError(
        "Pylinac is only supported on Python 3.7+. Please update your environment."
    )

# alphabetized modules
from .acr import ACRCT, ACRMRILarge
from .calibration import tg51, trs398
from .cheese import CIRS062M, TomoCheese

# import shortcuts
# core first
from .core import decorators, geometry, image, io, mask, profile, roi, utilities
from .core.profile import Centering
from .core.utilities import assign2machine, clear_data_files
from .ct import CatPhan503, CatPhan504, CatPhan600, CatPhan604
from .field_analysis import (
    Device,
    DeviceFieldAnalysis,
    Edge,
    FieldAnalysis,
    Interpolation,
    Normalization,
    Protocol,
)
from .log_analyzer import Dynalog, MachineLogs, TrajectoryLog, load_log
from .picketfence import PicketFence  # must be after log analyzer
from .planar_imaging import (
    PTWEPIDQC,
    SNCFSQA,
    SNCMV,
    SNCMV12510,
    DoselabMC2kV,
    DoselabMC2MV,
    DoselabRLf,
    ElektaLasVegas,
    IBAPrimusA,
    IMTLRad,
    IsoAlign,
    LasVegas,
    LeedsTOR,
    LeedsTORBlue,
    SNCkV,
    StandardImagingFC2,
    StandardImagingQC3,
    StandardImagingQCkV,
)
from .quart import HypersightQuartDVT, QuartDVT
from .starshot import Starshot
from .vmat import DRGS, DRMLC
from .winston_lutz import WinstonLutz, WinstonLutzMultiTargetMultiField
