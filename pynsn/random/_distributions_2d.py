from abc import ABCMeta, abstractmethod
from typing import Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray, ArrayLike

from . import _rng
from ..types import NoSolutionError

from ._distributions import ABCDistribution, round_samples


class MultiVarDistributionType(ABCDistribution):

    def __init__(self,
                 xy_minmax: Optional[ArrayLike] = None,
                 max_radius: Optional[float] = None):
        """
        xy_minmax: [xmin, ymin, xmax, ymax]
        """

        if xy_minmax is not None:
            xy_minmax = np.asarray(xy_minmax)
            if len(xy_minmax) != 4:
                raise ValueError("xy_minmax has to be an array of four values [xmin, ymin, xmax, ymax], "
                                 f"and not {xy_minmax}")

            # replace nan or None with np.inf/-np.inf
            for i, x in enumerate([-1*np.inf, -1*np.inf, np.inf, np.inf]):
                if xy_minmax[i] is None or np.isnan(xy_minmax[i]):
                    xy_minmax[i] = x

            if (xy_minmax[0] > xy_minmax[2] or xy_minmax[1] > xy_minmax[3]):
                raise TypeError(f"xy_minmax=[xmin, ymin, xmax, ymax]. xmin (or ymin) "
                                "is larger xmax (or ymax)")

        self.xy_minmax = xy_minmax
        self.max_radius = max_radius

    def to_dict(self):
        """Dict representation of the distribution"""
        d = super().to_dict()
        d.update({"max_radius": self.max_radius,
                  "xy_minmax": self.xy_minmax})
        return d


class Normal2D(MultiVarDistributionType):

    def __init__(self, mu: Tuple[float, float],
                 sigma: Tuple[float, float],
                 correlation: float = 0,
                 xy_minmax: Optional[ArrayLike] = None,
                 max_radius: Optional[float] = None):
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
        super().__init__(xy_minmax=xy_minmax, max_radius=max_radius)
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

            if self.xy_minmax is not None:
                idx = (draw[:, 0] < self.xy_minmax[0]) | \
                    (draw[:, 1] < self.xy_minmax[1]) | \
                    (draw[:, 0] > self.xy_minmax[2]) | \
                    (draw[:, 1] > self.xy_minmax[3])
                draw = np.delete(draw, idx, axis=0)

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
