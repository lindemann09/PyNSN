from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from copy import copy
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._shapes.colour import Colour
from .._shapes.shapes import PolygonShape
from ..exceptions import NoSolutionError
from . import _rng

ConstantLike = float | int | str | dict | PolygonShape | Colour
CategoricalLike = Sequence | NDArray


class AbstractDistribution(metaclass=ABCMeta):
    """Base class for all distribution"""

    @abstractmethod
    def to_dict(self) -> dict:
        """Dict representation of the distribution"""
        return {"type": type(self).__name__}

    @abstractmethod
    def sample(self, n: int) -> NDArray:
        """Random sample from the distribution

        Args:
            n: number of samples

        Returns:
            Numpy array of the sample

        """


class AbstractUnivarDistr(AbstractDistribution, metaclass=ABCMeta):
    pass


class AbstractContinuousDistr(AbstractUnivarDistr, metaclass=ABCMeta):
    """Univariate Continuous Distribution"""

    def __init__(self, minmax: ArrayLike | None):
        if minmax is None:
            self._minmax = np.array((None, None))
        else:
            # FIXME make properties with setter!
            self._minmax = np.asarray(minmax)
        if len(self._minmax) != 2:
            raise TypeError(f"min_max {minmax} has to be a tuple of two values")

    def to_dict(self) -> dict:
        """Dict representation of the distribution"""
        d = super().to_dict()
        d.update({"minmax": self._minmax.tolist()})
        return d

    @property
    def minmax(self) -> NDArray:
        return self._minmax


class Uniform(AbstractContinuousDistr):
    """ """

    def __init__(self, minmax: ArrayLike):
        """Uniform distribution defined by the number range, min_max=(min, max)

        Args:
            min_max : tuple (numeric, numeric)
                the range of the distribution
        """

        super().__init__(minmax)
        if self._minmax[0] > self._minmax[1]:
            raise TypeError(
                f"min_max {minmax} has to be a tuple of two values "
                "(a, b) with a <= b."
            )
        self._scale = self._minmax[1] - self._minmax[0]

    def sample(self, n: int) -> NDArray[np.float64]:
        dist = _rng.generator.random(size=n)
        return self._minmax[0] + dist * self._scale


class Triangle(AbstractContinuousDistr):
    """Triangle"""

    def __init__(self, mode: float, minmax: ArrayLike):
        super().__init__(minmax=minmax)
        self._mode = mode
        if (self._minmax[0] is not None and mode <= self._minmax[0]) or (
            self._minmax[1] is not None and mode >= self._minmax[1]
        ):
            raise ValueError(
                f"mode ({mode}) has to be inside the defined "
                f"min_max range ({self._minmax})"
            )

    def sample(self, n: int) -> NDArray[np.float64]:
        return _rng.generator.triangular(
            left=self._minmax[0], right=self._minmax[1], mode=self._mode, size=n
        )

    def to_dict(self) -> dict:
        d = super().to_dict()
        d.update({"mode": self._mode})
        return d

    @property
    def mode(self) -> float:
        return self._mode


class _AbstractDistrMuSigma(AbstractContinuousDistr, metaclass=ABCMeta):
    def __init__(self, mu: float, sigma: float, minmax: ArrayLike | None = None):
        super().__init__(minmax)
        self._mu = mu
        self._sigma = abs(sigma)
        if (self._minmax[0] is not None and mu <= self._minmax[0]) or (
            self._minmax[1] is not None and mu >= self._minmax[1]
        ):
            raise ValueError(
                f"mode ({mu}) has to be inside the defined "
                f"min_max range ({self._minmax})"
            )

    def to_dict(self) -> dict:
        d = super().to_dict()
        d.update({"mu": self._mu, "sigma": self._sigma})
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

    def sample(self, n: int) -> NDArray[np.float64]:
        rtn = np.array([])
        required = n
        while required > 0:
            draw = _rng.generator.normal(loc=self._mu, scale=self._sigma, size=required)
            if self._minmax[0] is not None:
                draw = np.delete(draw, draw < self._minmax[0])
            if self._minmax[1] is not None:
                draw = np.delete(draw, draw > self._minmax[1])
            if len(draw) > 0:  # type: ignore
                rtn = np.append(rtn, draw)
                required = n - len(rtn)
        return rtn


class Beta(_AbstractDistrMuSigma):
    def __init__(
        self,
        mu=None,
        sigma=None,
        alpha=None,
        beta=None,
        minmax: ArrayLike | None = None,
    ):
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
            raise TypeError("Either Mu & Sigma or Alpha & Beta have to specified.")
        super().__init__(mu=mu, sigma=sigma, minmax=minmax)

    def sample(self, n: int) -> NDArray[np.float64]:
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
        x = (m * (1 - m) / std**2) - 1
        return x * m, (1 - m) * x

    @property
    def alpha(self):
        return self.shape_parameter[0]

    @property
    def beta(self):
        return self.shape_parameter[1]

    @staticmethod
    def _calc_mu_sigma(
        alpha: float, beta: float, min_max: NDArray
    ) -> Tuple[float, float]:
        a = alpha
        b = beta
        r = min_max[1] - min_max[0]

        e = a / (a + b)
        mu = e * r + min_max[0]

        v = (a * b) / ((a + b) ** 2 * (a + b + 1))
        sigma = np.sqrt(v) * r
        return mu, sigma


class Categorical(AbstractUnivarDistr):
    """Categorical"""

    def __init__(
        self,
        levels: CategoricalLike,
        weights: ArrayLike | None = None,
        exact_weighting=False,
    ):
        """Distribution of category. Samples from discrete categories
        with optional weights for each category or category.
        """

        self._levels = np.asarray(copy(levels))
        self.exact_weighting = exact_weighting
        if weights is None:
            self._weights = np.empty(0)
        else:
            self._weights = np.asarray(weights)
            if len(self._levels) != len(self._weights):
                raise ValueError(
                    "Number weights does not match the number of category levels"
                )

    @property
    def levels(self) -> NDArray:
        return self._levels

    @property
    def weights(self) -> NDArray:
        return self._weights

    def sample(self, n: int) -> NDArray[np.float64]:
        if len(self._weights) == 0:
            p = np.array([1 / len(self._levels)] * len(self._levels))
        else:
            p = self._weights / np.sum(self._weights)

        if not self.exact_weighting:
            dist = _rng.generator.choice(a=self._levels, p=p, size=n)
        else:
            n_distr = n * p
            if np.any(np.round(n_distr) != n_distr):
                # problem: some n are floats
                try:
                    # greatest common denominator
                    gcd = np.gcd.reduce(self._weights)
                    info = "\nSample size has to be a multiple of {}.".format(
                        int(np.sum(self._weights / gcd))
                    )
                except:
                    info = ""

                raise NoSolutionError(
                    f"Can't find n={n} samples that"
                    + f" are exactly distributed as specified by the weights (p={p}). "
                    + info
                )

            dist = []
            for lev, n in zip(self._levels, n_distr):
                dist.extend([lev] * int(n))
            _rng.generator.shuffle(dist)

        return np.asarray(dist)

    def to_dict(self) -> dict:
        d = super().to_dict()
        d.update(
            {
                "levels": self._levels.tolist(),
                "weights": self._weights.tolist(),
                "exact": self.exact_weighting,
            }
        )
        return d


class Constant(AbstractUnivarDistr):
    def __init__(self, value: ConstantLike) -> None:
        """Helper class to "sample" constance.

        Looks like a PyNSNDistribution, but sample returns just the constant

        Parameter:
        ----------
        constant : numeric
        """

        self.value = value

    def sample(self, n: int) -> NDArray:
        return np.full(n, self.value)

    def to_dict(self) -> dict:
        return {"type": "Constant", "value": self.value}


DistributionLike = AbstractUnivarDistr | ConstantLike | CategoricalLike
