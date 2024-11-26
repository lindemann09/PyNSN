from abc import ABCMeta

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _rng
from ._distributions import AbstractDistribution


class Abstract2dDistr(AbstractDistribution, metaclass=ABCMeta):
    """Abstract class for Multivariate (2D) distributions"""

    def __init__(
        self,
        mu: ArrayLike,
        x_minmax: ArrayLike | None,
        y_minmax: ArrayLike | None,
        max_radius: float | None,
    ) -> None:
        """
        mu:
        x_range: array with one negative and one positive value
        y_range: relative bottom, top
        """

        self._mu = np.zeros(2)
        self._x_minmax = np.array((-1 * np.inf, np.inf))
        self._y_minmax = np.array((-1 * np.inf, np.inf))
        self._max_radius = np.inf
        # call constructors
        self.max_radius = max_radius
        self.set_limits(mu=mu, x_minmax=x_minmax, y_minmax=y_minmax)

    @property
    def mu(self) -> NDArray:
        return self._mu

    @property
    def x_min(self) -> ArrayLike:
        return self._x_minmax

    @property
    def y_min(self) -> ArrayLike:
        return self._y_minmax

    def set_limits(
        self, mu: ArrayLike, x_minmax: ArrayLike | None, y_minmax: ArrayLike | None
    ) -> None:
        self._mu = np.asarray(mu)
        if len(self._mu) != 2:
            raise TypeError("Mu has to be an array of two values.")

        if x_minmax is None:
            self._x_minmax = np.array((-1 * np.inf, np.inf))
        else:
            xr = np.asarray(x_minmax)
            if len(xr) != 2:
                raise TypeError("x_range has to be an array of two values")
            if xr[0] is None or np.isnan(xr[0]):
                xr[0] = -1 * np.inf
            if xr[1] is None or np.isnan(xr[1]):
                xr[1] = np.inf
            if not (xr[0] <= self._mu[0] <= xr[1]):
                raise ValueError(f"interval x_minmax {xr} must contain mu[0]")
            self._x_minmax = xr
        if y_minmax is None:
            self._y_minmax = np.array((-1 * np.inf, np.inf))
        else:
            yr = np.asarray(y_minmax)
            if len(yr) != 2:
                raise TypeError("y_range has to be an array of two values")
            if yr[0] is None or np.isnan(yr[0]):
                yr[0] = -1 * np.inf
            if yr[1] is None or np.isnan(yr[1]):
                yr[1] = np.inf
            if not (yr[0] <= self._mu[1] <= yr[1]):
                raise ValueError(f"interval y_minmax {yr} must contain mu[1]")
            self._y_minmax = yr

    @property
    def max_radius(self) -> float:
        return self._max_radius

    @max_radius.setter
    def max_radius(self, val: float | None):
        if val is None or np.isnan(val):
            self._max_radius = np.inf
        elif val <= 0:
            raise ValueError("max_radius has to be a positive float and not 0")
        else:
            self._max_radius = val

    def to_dict(self) -> dict:
        """Dict representation of the distribution"""
        d = super().to_dict()
        d.update(
            {
                "mu": self._mu.tolist(),
                "max_radius": self._max_radius,
                "x_minmax": self._x_minmax.tolist(),
                "y_minmax": self._y_minmax.tolist(),
            }
        )
        return d

    def _delete_outlier_cartesian(self, arr: NDArray):
        """helper: delete those 2D values that are outside a particular range
        (cartesian)
        """
        idx = (
            (arr[:, 0] < self._x_minmax[0])
            | (arr[:, 1] < self._y_minmax[0])
            | (arr[:, 0] > self._x_minmax[1])
            | (arr[:, 1] > self._y_minmax[1])
        )
        return np.delete(arr, idx, axis=0)

    def _delete_outlier_radial(self, arr: NDArray):
        """helper: delete those 2D values that are outside the radial area
        (mu, max_radius)
        """

        if self._max_radius < np.inf:
            # remove to large radii
            d = arr - self._mu
            r = np.hypot(d[:, 0], d[:, 1])
            return np.delete(arr, r > self._max_radius, axis=0)
        else:
            return arr


class Uniform2D(Abstract2dDistr):
    def __init__(
        self,
        x_minmax: ArrayLike | None = None,
        y_minmax: ArrayLike | None = None,
        radial_mu: ArrayLike | None = None,
        max_radius: float | None = None,
    ):
        """Two dimensional uniform distribution with cut-off either cardinal or
        radial

        Args:
            x_minmax, y_minmax: for cardinal limits
            radial_mu, max_radius: for radial limits
        """
        cardinal = x_minmax is not None and y_minmax is not None
        radial = radial_mu is not None and max_radius is not None
        if cardinal == radial:  # cardinal & radial or both false
            raise ValueError("Define either x_min & ymin OR radial_mu & radial_radius")
        elif cardinal:
            x_minmax = np.asarray(x_minmax)
            y_minmax = np.asarray(y_minmax)
            mu = np.array([np.mean(x_minmax), np.mean(y_minmax)])
            max_radius = None
        else:
            # radial
            mu = np.asarray(radial_mu)
            r = np.array((-1 * max_radius, max_radius))  # type: ignore
            x_minmax = mu[0] + r
            y_minmax = mu[1] + r

        super().__init__(
            mu=mu, x_minmax=x_minmax, y_minmax=y_minmax, max_radius=max_radius
        )

        self._xy_scale = (
            self._x_minmax[1] - self._x_minmax[0],
            self._y_minmax[1] - self._y_minmax[0],
        )

    def sample(self, n: int) -> NDArray[np.float64]:
        rtn = np.empty((0, 2), dtype=float)
        required = n
        while required > 0:
            draw = _rng.generator.random(size=(required, 2))
            draw[:, 0] = self._x_minmax[0] + draw[:, 0] * self._xy_scale[0]
            draw[:, 1] = self._y_minmax[0] + draw[:, 1] * self._xy_scale[1]
            draw = self._delete_outlier_radial(draw)
            if len(draw) > 0:
                # append
                rtn = np.append(rtn, draw, axis=0)
                required = n - len(rtn)

        return rtn


class Normal2D(Abstract2dDistr):
    def __init__(
        self,
        mu: ArrayLike,
        sigma: ArrayLike,
        correlation: float = 0,
        x_minmax: ArrayLike | None = None,
        y_minmax: ArrayLike | None = None,
        max_radius: float | None = None,
    ):
        """Two dimensional normal distribution with optional cut-off (radial or
        cartesian)

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
        super().__init__(
            mu=mu, x_minmax=x_minmax, y_minmax=y_minmax, max_radius=max_radius
        )
        self.correlation = correlation
        self.sigma = np.abs(sigma)

        if len(self.sigma) != 2:
            raise TypeError("sigma must be an array of two values")
        if correlation < -1 or correlation > 1:
            raise ValueError("Correlations has to be between -1 and 1")

    def varcov(self):
        """Variance covariance matrix"""
        v = self.sigma**2
        cov = np.eye(2) * v
        cov[0, 1] = self.correlation * np.sqrt(v[0] * v[1])
        cov[1, 0] = cov[0, 1]
        return cov

    def sample(self, n: int) -> NDArray[np.float64]:
        rtn = None
        required = n
        while required > 0:
            draw = _rng.generator.multivariate_normal(
                mean=self._mu, cov=self.varcov(), size=required
            )
            draw = self._delete_outlier_radial(self._delete_outlier_cartesian(draw))
            if len(draw) > 0:
                # append
                if rtn is None:
                    rtn = draw
                else:
                    rtn = np.append(rtn, draw, axis=0)
                required = n - len(rtn)

        if rtn is None:
            return np.empty(0, dtype=float)
        else:
            return rtn

    def to_dict(self):
        d = super().to_dict()
        d.update({"sigma": self.sigma.tolist(), "correlation": self.correlation})
        return d
