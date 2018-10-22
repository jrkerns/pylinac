from typing import Union

import numpy as np


NumberLike = Union[int, np.float, float]

ArrayLike = Union[list, tuple, np.ndarray]

NumberOrArray = Union[NumberLike, ArrayLike]