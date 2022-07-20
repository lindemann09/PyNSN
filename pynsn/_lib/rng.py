"""NOTE:
import always ethe ntire module and call `rng.generate` to ensure always the
access of the newly initialized random generator if  `init_random_generator` has
been called
"""

import numpy as np
from .lib_typing import Union, ArrayLike
generator = np.random.default_rng()


def init_random_generator(seed: Union[int, ArrayLike, None] = None):
    """Init random generator and set random seed (optional)

    Parameters
    ----------
    seed: seed value
        must be none, int or array_like[ints]

    Notes
    -----
    see documentation of `numpy.random.default_rng()` of Python standard library
    """
    global generator
    generator = np.random.default_rng(seed=seed)
    if seed is not None:
        print("PyNSN seed: {}".format(seed))
