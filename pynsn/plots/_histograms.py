from matplotlib.pyplot import hist, hist2d

from ..rnd._distributions import AbstractDistribution


def distribution_samples(distr: AbstractDistribution,
                         n: int = 100000):
    """Creating a visualization of the distribution with ``matplotlib.pyplot``

    Parameters
    ----------
    distr : Abstract2dDistr or AbstractUnivarDistr
        distribution object
    n : int, optional
         number of sample of samples, by default 100000

    Returns
    -------
    ``matplotlib.pyplot.hist()`` or ``matplotlib.pyplot.hist2d()``

    Notes
    -----
    call plt.show() to display the figure

    """

    samples = distr.sample(n=n)
    if samples.ndim == 2:
        return hist2d(samples[:, 0], samples[:, 1], bins=(100, 100))[2]
    else:
        return hist(samples, bins=100)[2]
