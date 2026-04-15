"""Typing utilities for Pylinac."""

from __future__ import annotations

import numpy as np

ArrayLike = list | tuple | np.ndarray

NumberOrArray = float | ArrayLike
