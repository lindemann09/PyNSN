from abc import ABCMeta, abstractmethod
from typing import Sequence, Union

import numpy as np
from numpy.typing import NDArray, ArrayLike

from . import _rng
from ..types import NoSolutionError

from ._distributions import ABCDistribution, round_samples


class MultiVarDistributionType(ABCDistribution):

    def __init__(self, x_minmax, y_minmax):

        self.x_minmax = np.asarray(x_minmax)
        self.y_minmax = np.asarray(y_minmax)
        if len(self.x_minmax) != 2 or \
                (None not in self.x_minmax and self.x_minmax[0] > self.x_minmax[1]):
            raise TypeError(f"min_max {x_minmax} has to be a tuple of two values "
                            "(a, b) with a <= b.")
        if len(self.y_minmax) != 2 or \
                (None not in self.y_minmax and self.y_minmax[0] > self.y_minmax[1]):
            raise TypeError(f"min_max {y_minmax} has to be a tuple of two values "
                            "(a, b) with a <= b.")

    def to_dict(self):
        """Dict representation of the distribution"""
        d = super().to_dict()
        d.update({"x_min_max": self.x_minmax,
                  "y_min_max": self.y_minmax})
        return d


class Normal2D(MultiVarDistributionType):

    def __init__(self, mu, sigma, correlation, max_radius=None):
        """Two dimensional normal distribution with optional radial cut-off

        Resulting multivariate distribution has the defined means, standard
        deviations and correlation (approx.) between the two variables

        Args:
            mu:  tuple of to numbers (numeric, numeric)
            sigma: tuple of to numbers  (numeric, numeric)
            correlation: numeric
                the correlation between a and b
            max_radius: numerical (optional)
                cut-off criterion: maximal distance from the distribution center (mu)
        """
        super().__init__(min_max=np.array((None, None)))

        try:
            lmu = len(mu)
            lsigma = len(sigma)
        except TypeError:
            lmu = None
            lsigma = None
        if lmu != 2 or lsigma != 2:
            raise ValueError("Mu and sigma has to be a tuple of two values.")
        if correlation < -1 or correlation > 1:
            raise ValueError("Correlations has to be between -1 and 1")

        self.mu = np.array(mu)
        self.sigma = np.abs(sigma)
        self.correlation = correlation
        self.max_radius = max_radius

    def varcov(self):
        """Variance covariance matrix"""
        v = self.sigma ** 2
        cov = np.eye(2) * v
        cov[0, 1] = self.correlation * np.sqrt(v[0] * v[1])
        cov[1, 0] = cov[0, 1]
        return cov

    def sample(self, n, round_to_decimals=None):
        rtn = None
        required = n
        while required > 0:
            draw = _rng.generator.multivariate_normal(
                mean=self.mu, cov=self.varcov(), size=required)
            if self.max_radius is not None:
                # remove to large radii
                r = np.hypot(draw[:, 0], draw[:, 1])
                draw = np.delete(draw, r > self.max_radius, axis=0)

            if len(draw) > 0:
                # append
                if rtn is None:
                    rtn = draw
                else:
                    rtn = np.append(rtn, draw, axis=0)
                required = n - len(rtn)

        if rtn is None:
            return None
        else:
            return round_samples(rtn, round_to_decimals)

    def to_dict(self):
        d = super().to_dict()
        d.update({"mu": self.mu.tolist(),
                  "sigma": self.sigma.tolist(),
                  "correlation": self.correlation,
                  "max_radius": self.max_radius})
        return d
