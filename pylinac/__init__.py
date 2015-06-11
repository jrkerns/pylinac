
__version__ = '0.6.0.5'
__version_info__ = (0, 6, 0, 5)

import sys
if sys.version_info[0] < 3:
    raise ValueError("Pylinac is only supported on Python 3.x. It seems you are using Python 2; please use a different interpreter.")

DEBUG = False
MEMORY_PROFILE = False

if DEBUG:
    import warnings
    warnings.warn("Pylinac is in Debug mode.")