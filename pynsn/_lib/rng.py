import numpy as np

generator = np.random.default_rng()


def init_random_generator(seed=None):
    """Init random generator and set random seed

    Parameters
    ----------
    a : seed value
        must be  int, array_like[ints], SeedSequence, BitGenerator, Generator

    Notes
    -----
    see documentation of `numpy.random.default_rng()` of Python standard library
    """
    global generator
    generator = np.random.default_rng(seed=seed)
    if seed is not None:
        print("seed: {}".format(seed))
