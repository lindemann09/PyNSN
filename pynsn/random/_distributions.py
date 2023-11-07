from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from typing import Dict, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..types import NoSolutionError
from . import _rng


def round_samples(samples: NDArray,
                  round_to_decimals: Union[None, int]) -> NDArray:
    if round_to_decimals is not None:
        arr = np.round(samples, decimals=round_to_decimals)
        if round_to_decimals == 0:
            return arr.astype(np.int32)
        else:
            return arr
    else:
        return np.asarray(samples)


class ABCDistribution(metaclass=ABCMeta):

    @abstractmethod
    def to_dict(self) -> Dict:
        """Dict representation of the distribution"""
        return {"distribution": type(self).__name__}

    @abstractmethod
    def sample(self, n, round_to_decimals=False) -> NDArray:
        """Random sample from the distribution

        Args:
            n: number of samples
            round_to_decimals: Set to round samples. If 0 a array of integer
                will be return

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
            return hist2d(samples[:, 0], samples[:, 1], bins=100)[2]
        else:
            return hist(samples, bins=100)[2]


class UnivariateDistributionType(ABCDistribution):
    """Univariate Distribution
    """

    def __init__(self, minmax: Optional[ArrayLike]):
        if minmax is None:
            self.minmax = np.array((None, None))
        else:
            self.minmax = np.asarray(minmax)
        if len(self.minmax) != 2:
            raise TypeError(
                f"min_max {minmax} has to be a tuple of two values")

    def to_dict(self):
        """Dict representation of the distribution"""
        d = super().to_dict()
        d.update({"min_max": self.minmax})
        return d


class Uniform(UnivariateDistributionType):
    """
    """

    def __init__(self, minmax: ArrayLike):
        """Uniform distribution defined by the number range, min_max=(min, max)

        Args:
            min_max : tuple (numeric, numeric)
                the range of the distribution
        """

        super().__init__(minmax)
        if self.minmax[0] is None or self.minmax[1] is None or \
                self.minmax[0] > self.minmax[1]:
            raise TypeError(f"min_max {minmax} has to be a tuple of two values "
                            "(a, b) with a <= b.")

    def sample(self, n, round_to_decimals=None):
        dist = _rng.generator.random(size=n)
        rtn = self.minmax[0] + dist * float(self.minmax[1] - self.minmax[0])
        return round_samples(rtn, round_to_decimals)


class Levels(UnivariateDistributionType):
    """Levels
    """

    def __init__(self, levels: Sequence,
                 weights: Optional[ArrayLike] = None,
                 exact_weighting=False):
        """Distribution of level. Samples from a population discrete categories
         with optional weights for each level or category.
        """
        super().__init__(minmax=None)
        self.levels = levels
        self.exact_weighting = exact_weighting
        if weights is None:
            self.weights = None
        else:
            self.weights = np.asarray(weights)
            if len(levels) != len(self.weights):
                raise ValueError(
                    "Number weights does not match the number of levels")

    def sample(self, n, round_to_decimals=None):
        if self.weights is None:
            p = np.array([1 / len(self.levels)] * len(self.levels))
        else:
            p = self.weights / np.sum(self.weights)

        if not self.exact_weighting:
            dist = _rng.generator.choice(a=self.levels, p=p, size=n)
        else:
            n_distr = n * p
            if np.any(np.round(n_distr) != n_distr):
                # problem: some n are floats
                try:
                    # greatest common denominator
                    gcd = np.gcd.reduce(self.weights)  # type: ignore
                    info = "\nSample size has to be a multiple of {}.".format(
                        int(np.sum(self.weights / gcd)))
                except:
                    info = ""

                raise NoSolutionError(f"Can't find n={n} samples that" +
                                      f" are exactly distributed as specified by the weights (p={p}). " +
                                      info)

            dist = []
            for lev, n in zip(self.levels, n_distr):
                dist.extend([lev] * int(n))
            _rng.generator.shuffle(dist)

        return round_samples(np.asarray(dist), round_to_decimals)

    def to_dict(self):
        d = super().to_dict()
        d.update({"population": self.levels,
                  "weights": self.weights,
                  "exact_weighting": self.exact_weighting})
        return d


class Triangle(UnivariateDistributionType):
    """Triangle
    """

    def __init__(self, mode, minmax: ArrayLike):
        super().__init__(minmax=minmax)
        self.mode = mode
        if (self.minmax[0] is not None and mode <= self.minmax[0]) or \
                (self.minmax[1] is not None and mode >= self.minmax[1]):
            raise ValueError(f"mode ({mode}) has to be inside the defined "
                             f"min_max range ({self.minmax})")

    def sample(self, n, round_to_decimals=None):
        dist = _rng.generator.triangular(left=self.minmax[0], right=self.minmax[1],
                                         mode=self.mode, size=n)
        return round_samples(dist, round_to_decimals)

    def to_dict(self):
        d = super().to_dict()
        d.update({"mode": self.mode})
        return d


class _DistributionMuSigma(UnivariateDistributionType):

    def __init__(self, mu, sigma, minmax: Optional[ArrayLike] = None):
        super().__init__(minmax)
        self.mu = mu
        self.sigma = abs(sigma)
        if (self.minmax[0] is not None and mu <= self.minmax[0]) or \
                (self.minmax[1] is not None and mu >= self.minmax[1]):
            raise ValueError(f"mode ({mu}) has to be inside the defined "
                             f"min_max range ({self.minmax})")

    def to_dict(self):
        d = super().to_dict()
        d.update({"mu": self.mu,
                  "sigma": self.sigma})
        return d


class Normal(_DistributionMuSigma):
    """Normal distribution with optional cut-off of minimum and maximum

    Resulting distribution has the defined mean and std, only if
    cutoffs values are symmetric.

    Args:
        mu: numeric
        sigma: numeric
        min_max : tuple (numeric, numeric) or None
            the range of the distribution
    """

    def sample(self, n, round_to_decimals=None):
        rtn = np.array([])
        required = n
        while required > 0:
            draw = _rng.generator.normal(
                loc=self.mu, scale=self.sigma, size=required)
            if self.minmax[0] is not None:
                draw = np.delete(draw, draw < self.minmax[0])
            if self.minmax[1] is not None:
                draw = np.delete(draw, draw > self.minmax[1])
            if len(draw) > 0:  # type: ignore
                rtn = np.append(rtn, draw)
                required = n - len(rtn)
        return round_samples(rtn, round_to_decimals)


class Beta(_DistributionMuSigma):

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

    def sample(self, n, round_to_decimals=None):
        if self.sigma is None or self.sigma == 0:
            return np.array([self.mu] * n)

        alpha, beta = self.shape_parameter
        dist = _rng.generator.beta(a=alpha, b=beta, size=n)
        dist = (dist - np.mean(dist)) / np.std(dist)  # z values
        rtn = dist * self.sigma + self.mu
        return round_samples(rtn, round_to_decimals)

    @property
    def shape_parameter(self):
        """Alpha (p) & beta (q) parameter for the beta distribution
        http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm

        Returns
        -------
        parameter: tuple
            shape parameter (alpha, beta) of the distribution

        """
        r = float(self.minmax[1] - self.minmax[0])
        m = (self.mu - self.minmax[0]) / r  # mean
        std = self.sigma / r
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


class _Constant(UnivariateDistributionType):

    def __init__(self, value: float) -> None:
        """Helper class to "sample" constance.

        Looks like a PyNSNDistribution, but sample returns just the constant

        Parameter:
        ----------
        constant : numeric
        """

        super().__init__(minmax=np.array((value, value)))

    def sample(self, n, round_to_decimals=None) -> NDArray:
        return round_samples([self.minmax[0]] * n, round_to_decimals)

    def to_dict(self) -> dict:
        return {"distribution": "Constant",
                "value": self.minmax[0]}
