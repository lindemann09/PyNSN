from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from copy import copy
from typing import Any, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..errors import NoSolutionError
from . import _rng


class AbstractDistribution(metaclass=ABCMeta):
    """Base class for all distribution"""

    @abstractmethod
    def todict(self) -> dict:
        """Dict representation of the distribution"""
        return {"type": type(self).__name__}

    @abstractmethod
    def sample(self, n: int) -> NDArray[np.float_]:
        """Random sample from the distribution

        Args:
            n: number of samples

        Returns:
            Numpy array of the sample

        """

    def pyplot_samples(self, n=100000):
        """Creating a visualization of the distribution with ``matplotlib.pyplot``

        Args:
            n: number of sample (optional)

        Returns:
            ``pyplot.figure()``

        Raises:
            ImportError: if ``matplotlib`` is not installed
        """
        try:
            from matplotlib.pyplot import hist, hist2d
        except ImportError as err:
            raise ImportError(
                "To use pyplot_samples, please install matplotlib.") from err

        samples = self.sample(n=n)

        if samples.ndim == 2:
            return hist2d(samples[:, 0], samples[:, 1], bins=(100, 100))[2]
        else:
            return hist(samples, bins=100)[2]


class AbstractUnivarDistr(AbstractDistribution, metaclass=ABCMeta):
    """Univariate Distribution
    """

    def __init__(self, minmax: Optional[ArrayLike]):
        if minmax is None:
            self._minmax = np.array((None, None))
        else:
            # FIXME make properties with setter!
            self._minmax = np.asarray(minmax)
        if len(self._minmax) != 2:
            raise TypeError(
                f"min_max {minmax} has to be a tuple of two values")

    def todict(self) -> dict:
        """Dict representation of the distribution"""
        d = super().todict()
        d.update({"minmax": self._minmax.tolist()})
        return d

    @property
    def minmax(self) -> NDArray:
        return self._minmax


class Uniform(AbstractUnivarDistr):
    """
    """

    def __init__(self, minmax: ArrayLike):
        """Uniform distribution defined by the number range, min_max=(min, max)

        Args:
            min_max : tuple (numeric, numeric)
                the range of the distribution
        """

        super().__init__(minmax)
        if self._minmax[0] > self._minmax[1]:
            raise TypeError(f"min_max {minmax} has to be a tuple of two values "
                            "(a, b) with a <= b.")
        self._scale = self._minmax[1] - self._minmax[0]

    def sample(self, n: int) -> NDArray[np.float_]:
        dist = _rng.generator.random(size=n)
        return self._minmax[0] + dist * self._scale


class Triangle(AbstractUnivarDistr):
    """Triangle
    """

    def __init__(self, mode: float, minmax: ArrayLike):
        super().__init__(minmax=minmax)
        self._mode = mode
        if (self._minmax[0] is not None and mode <= self._minmax[0]) or \
                (self._minmax[1] is not None and mode >= self._minmax[1]):
            raise ValueError(f"mode ({mode}) has to be inside the defined "
                             f"min_max range ({self._minmax})")

    def sample(self, n: int) -> NDArray[np.float_]:
        return _rng.generator.triangular(left=self._minmax[0], right=self._minmax[1],
                                         mode=self._mode, size=n)

    def todict(self) -> dict:
        d = super().todict()
        d.update({"mode": self._mode})
        return d

    @property
    def mode(self) -> float:
        return self._mode


class _AbstractDistrMuSigma(AbstractUnivarDistr, metaclass=ABCMeta):

    def __init__(self, mu: float, sigma: float, minmax: Optional[ArrayLike] = None):
        super().__init__(minmax)
        self._mu = mu
        self._sigma = abs(sigma)
        if (self._minmax[0] is not None and mu <= self._minmax[0]) or \
                (self._minmax[1] is not None and mu >= self._minmax[1]):
            raise ValueError(f"mode ({mu}) has to be inside the defined "
                             f"min_max range ({self._minmax})")

    def todict(self) -> dict:
        d = super().todict()
        d.update({"mu": self._mu,
                  "sigma": self._sigma})
        return d

    @property
    def mu(self) -> float:
        return self._mu

    @property
    def sigma(self) -> float:
        return self._sigma


class Normal(_AbstractDistrMuSigma):
    """Normal distribution with optional cut-off of minimum and maximum

    Resulting distribution has the defined mean and std, only if
    cutoffs values are symmetric.

    Args:
        mu: numeric
        sigma: numeric
        min_max : tuple (numeric, numeric) or None
            the range of the distribution
    """

    def sample(self, n: int) -> NDArray[np.float_]:
        rtn = np.array([])
        required = n
        while required > 0:
            draw = _rng.generator.normal(
                loc=self._mu, scale=self._sigma, size=required)
            if self._minmax[0] is not None:
                draw = np.delete(draw, draw < self._minmax[0])
            if self._minmax[1] is not None:
                draw = np.delete(draw, draw > self._minmax[1])
            if len(draw) > 0:  # type: ignore
                rtn = np.append(rtn, draw)
                required = n - len(rtn)
        return rtn


class Beta(_AbstractDistrMuSigma):

    def __init__(self, mu=None, sigma=None, alpha=None, beta=None,
                 minmax: Optional[ArrayLike] = None):
        """Beta distribution defined by the number range, min_max=(min, max),
         the mean and the standard deviation (std)

        Resulting distribution has the defined mean and std

        Args:
            mu: numeric
            sigma: numeric
            min_max : tuple (numeric, numeric)
                the range of the distribution

        Note:
            Depending on the position of the mean in number range the
            distribution is left or right skewed.

        See Also:
            for calculated shape parameters [alpha, beta] see property
            `shape_parameter_beta`
        """
        if alpha is not None and beta is not None and (mu, sigma) == (None, None):
            minmax = np.asarray(minmax)
            mu, sigma = Beta._calc_mu_sigma(alpha, beta, minmax)
        elif mu is None or sigma is None or alpha is not None or beta is not None:
            raise TypeError(
                "Either Mu & Sigma or Alpha & Beta have to specified.")
        super().__init__(mu=mu, sigma=sigma, minmax=minmax)

    def sample(self, n: int) -> NDArray[np.float_]:
        if self._sigma is None or self._sigma == 0:
            return np.array([self._mu] * n)

        alpha, beta = self.shape_parameter
        dist = _rng.generator.beta(a=alpha, b=beta, size=n)
        dist = (dist - np.mean(dist)) / np.std(dist)  # z values
        rtn = dist * self._sigma + self._mu
        return rtn

    @property
    def shape_parameter(self):
        """Alpha (p) & beta (q) parameter for the beta distribution
        http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm

        Returns
        -------
        parameter: tuple
            shape parameter (alpha, beta) of the distribution

        """
        r = float(self._minmax[1] - self._minmax[0])
        m = (self._mu - self._minmax[0]) / r  # mean
        std = self._sigma / r
        x = (m * (1 - m) / std ** 2) - 1
        return x * m, (1 - m) * x

    @property
    def alpha(self):
        return self.shape_parameter[0]

    @property
    def beta(self):
        return self.shape_parameter[1]

    @staticmethod
    def _calc_mu_sigma(alpha: float, beta: float, min_max: NDArray) -> Tuple[float, float]:
        a = alpha
        b = beta
        r = min_max[1] - min_max[0]

        e = a / (a + b)
        mu = e * r + min_max[0]

        v = (a * b) / ((a + b) ** 2 * (a + b + 1))
        sigma = np.sqrt(v) * r
        return mu, sigma


class Categorical(AbstractDistribution):
    """Categorical
    """

    def __init__(self,
                 categories: Union[Sequence, NDArray],
                 weights: Optional[ArrayLike] = None,
                 exact_weighting=False):
        """Distribution of category. Samples from a population discrete categories
         with optional weights for each category or category.
        """

        self._categories = np.asarray(copy(categories))
        self.exact_weighting = exact_weighting
        if weights is None:
            self._weights = np.empty(0)
        else:
            self._weights = np.asarray(weights)
            if len(self._categories) != len(self._weights):
                raise ValueError(
                    "Number weights does not match the number of categories")

    @property
    def categories(self) -> NDArray:
        return self._categories

    @property
    def weights(self) -> NDArray:
        return self._weights

    def sample(self, n: int) -> NDArray[np.float_]:
        if len(self._weights) == 0:
            p = np.array([1 / len(self._categories)] * len(self._categories))
        else:
            p = self._weights / np.sum(self._weights)

        if not self.exact_weighting:
            dist = _rng.generator.choice(a=self._categories, p=p, size=n)
        else:
            n_distr = n * p
            if np.any(np.round(n_distr) != n_distr):
                # problem: some n are floats
                try:
                    # greatest common denominator
                    gcd = np.gcd.reduce(self._weights)
                    info = "\nSample size has to be a multiple of {}.".format(
                        int(np.sum(self._weights / gcd)))
                except:
                    info = ""

                raise NoSolutionError(f"Can't find n={n} samples that" +
                                      f" are exactly distributed as specified by the weights (p={p}). " +
                                      info)

            dist = []
            for lev, n in zip(self._categories, n_distr):
                dist.extend([lev] * int(n))
            _rng.generator.shuffle(dist)

        return np.asarray(dist)

    def todict(self) -> dict:
        d = super().todict()
        d.update({"population": self._categories.tolist(),
                  "weights": self._weights.tolist(),
                  "exact_weighting": self.exact_weighting})
        return d


class Constant(AbstractDistribution):

    def __init__(self, value: Any) -> None:
        """Helper class to "sample" constance.

        Looks like a PyNSNDistribution, but sample returns just the constant

        Parameter:
        ----------
        constant : numeric
        """

        self.value = value

    def sample(self, n: int) -> NDArray[np.float_]:
        return np.full(self.value,  n)

    def todict(self) -> dict:
        return {"type": "Constant",
                "value": self.value}
