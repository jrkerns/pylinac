
__version__ = '0.4.1'
__version_info__ = (0, 4, 1)

import sys
if sys.version_info[0] < 3:
    raise ValueError("Pylinac is only supported on Python 3.x. It seems you are using Python 2; please use a different interpreter.")

try:
    import pandas
except ImportError:
    import warnings
    warnings.warn("Future versions of pylinac will use the pandas package, but it is not installed", FutureWarning)

DEBUG = False
MEMORY_PROFILE = False

if DEBUG:
    import warnings
    warnings.warn("Pylinac is in Debug mode.")