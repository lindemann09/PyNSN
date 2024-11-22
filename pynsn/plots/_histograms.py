from typing import Tuple
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from ..rnd._distributions import AbstractDistribution


def distribution_samples(distr: AbstractDistribution,
                         n: int = 100000,
                         ax: Axes | None = None) -> Tuple[Figure | None, Axes]:
    """Creating a visualization of the distribution with ``matplotlib.pyplot``

    Parameters
    ----------
    distr : Abstract2dDistr or AbstractUnivarDistr
        distribution object
    n : int, optional
        number of sample of samples, by default 100000
    ax : `matplotlib.axes.Axes` or None, default None

    Returns
    -------
    Figure, Axes

    Notes
    -----
    call ``plt.show()`` to display the figures

    """

    if ax is None:
        ax = plt.figure().subplots()

    samples = distr.sample(n=n)
    if samples.ndim == 2:
        ax.hist2d(samples[:, 0], samples[:, 1], bins=(100, 100))
    else:
        ax.hist(samples, bins=100)

    return ax.get_figure(), ax
