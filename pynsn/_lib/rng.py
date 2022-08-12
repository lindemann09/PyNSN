"""NOTE:
import always entire module and call `rng.generator` to ensure always the
access of the newly initialized random generator if  `init_random_generator` has
been called
"""

from typing import Sequence, Union

import numpy as np
from numpy.typing import NDArray

generator = np.random.default_rng()


def init_random_generator(seed: Union[int, NDArray, Sequence, None] = None):
    """Init random generator and set random seed (optional)

    Parameters
    ----------
    seed: seed value
        must be none, int or array_like[ints]

    Notes
    -----
    see documentation of `numpy.random.default_rng()` of Python standard library
    """

    # pylint:disable = W0603
    global generator
    generator = np.random.default_rng(seed=seed)
    if seed is not None:
        print(f"PyNSN seed: {seed}")
