import warnings

from ..metrics.features import *  # noqa:F403
from ..metrics.image import *  # noqa:F403
from ..metrics.utils import *  # noqa:F403

warnings.warn(
    "This module has been moved to pylinac.metrics. Please import from there in the future.",
    DeprecationWarning,
)
