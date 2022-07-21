import numpy as _np
import numpy as np

from .._lib.misc import numpy_round2 as _numpy_round2
from .._lib import rng as _rng
from .._lib.exception import NoSolutionError as _NoSolutionError


def _round_samples(samples, round_to_decimals):
    if round_to_decimals is not None:
        return _numpy_round2(samples, decimals=round_to_decimals)
    else:
        return _np.array(samples)


class PyNSNDistribution(object):

    def __init__(self, min_max):
        if not isinstance(min_max, (list, tuple)) or len(min_max) != 2 or \
                (None not in min_max and min_max[0] > min_max[1]):
            raise TypeError("min_max {} ".format(min_max) + \
                            "has to be a tuple of two values (a, b) with a <= b.")
        self.min_max = min_max

    def as_dict(self):
        return {"distribution": type(self).__name__,
                "min_max": self.min_max}

    def sample(self, n, round_to_decimals=False):
        return NotImplementedError()

    def pyplot_samples(self, n=100000):
        try:
            from matplotlib.pyplot import hist, hist2d
        except:
            raise ImportError("To use pyplot_samples, please install matplotlib.")

        samples = self.sample(n=n)
        if samples.ndim == 2:
            return hist2d(samples[:,0], samples[:,1], bins=100)[2]
        else:
            return hist(samples, bins=100)[2]


class Uniform(PyNSNDistribution):
    """
    """

    def __init__(self, min_max):
        """Uniform distribution defined by the number range, min_max=(min, max)

        Parameter:
        ----------
        min_max : tuple (numeric, numeric)
            the range of the distribution
        """

        super().__init__(min_max)

    def sample(self, n, round_to_decimals=None):
        dist = _rng.generator.random(size=n)
        rtn = self.min_max[0] + dist * float(self.min_max[1] - self.min_max[0])
        return _round_samples(rtn, round_to_decimals)


class Levels(PyNSNDistribution):
    """
    """

    def __init__(self, levels, weights=None, exact_weighting=False):
        """Distribution of level. Samples from a population discrete categories
         with optional weights for each level or category.

        Parameter:
        ----------
        """
        super().__init__(min_max=(None, None))
        self.levels = levels
        self.exact_weighting = exact_weighting

        if weights is not None and len(levels) != len(weights):
            raise ValueError("Number weights does not match the number of levels")
        self.weights = weights

    def sample(self, n, round_to_decimals=None):
        if self.weights is None:
            p = _np.array([1 / len(self.levels)] * len(self.levels))
        else:
            p = _np.array(self.weights)
            p = p / _np.sum(p)

        if not self.exact_weighting:
            dist = _rng.generator.choice(a=self.levels, p=p, size=n)
        else:
            n_distr = n * p
            if not _np.alltrue(_np.round(n_distr) == n_distr):
                # problem: some n are floats
                try:
                    gcd = _np.gcd.reduce(self.weights)  # greatest common denominator
                    info = "\nSample size has to be a multiple of {}.".format(
                        int(_np.sum(self.weights / gcd)))
                except:
                    info = ""
                raise _NoSolutionError(f"Can't find n={n} samples that" +
                      f" are exactly distributed as specified by the weights (p={p}). " +
                      info)

            dist = []
            for lev, n in zip(self.levels, n_distr):
                dist.extend([lev] * int(n))
            _rng.generator.shuffle(dist)

        return _round_samples(dist, round_to_decimals)

    def as_dict(self):
        d = super().as_dict()
        d.update({"population": self.levels,
                  "weights": self.weights,
                  "exact_weighting": self.exact_weighting})
        return d


class Triangle(PyNSNDistribution):
    """
    """

    def __init__(self, mode, min_max):
        super().__init__(min_max=min_max)
        self.mode = mode
        if (min_max[0] is not None and mode <= min_max[0]) or \
                (min_max[1] is not None and mode >= min_max[1]):
            txt = "mode ({}) has to be inside the defined min_max range ({})".format(
                mode, min_max)
            raise ValueError(txt)

    def sample(self, n, round_to_decimals=None):
        dist = _rng.generator.triangular(left=self.min_max[0], right=self.min_max[1],
                                         mode=self.mode, size=n)
        return _round_samples(dist, round_to_decimals)

    def as_dict(self):
        d = super().as_dict()
        d.update({"mode": self.mode})
        return d


class _PyNSNDistributionMuSigma(PyNSNDistribution):

    def __init__(self, mu, sigma, min_max):
        super().__init__(min_max)
        self.mu = mu
        self.sigma = abs(sigma)
        if (min_max[0] is not None and mu <= min_max[0]) or \
                (min_max[1] is not None and mu >= min_max[1]):
            txt = "mean ({}) has to be inside the defined min_max range ({})".format(
                mu, min_max)
            raise ValueError(txt)

    def as_dict(self):
        d = super().as_dict()
        d.update({"mu": self.mu,
                  "sigma": self.sigma})
        return d


class Normal(_PyNSNDistributionMuSigma):

    def __init__(self, mu, sigma, min_max=None):
        """Normal distribution with optional cut-off of minimum and maximum

        Resulting distribution has the defined mean and std, only if
        cutoffs values are symmetric.

        Parameter:
        ----------
        mu: numeric
        sigma: numeric
        min_max : tuple (numeric, numeric) or None
            the range of the distribution
        """
        if min_max is None:
            min_max = (None, None)
        super().__init__(mu, sigma, min_max)

    def sample(self, n, round_to_decimals=None):
        rtn = _np.array([])
        required = n
        while required > 0:
            draw = _rng.generator.normal(loc=self.mu, scale=self.sigma, size=required)
            if self.min_max[0] is not None:
                draw = _np.delete(draw, draw < self.min_max[0])
            if self.min_max[1] is not None:
                draw = _np.delete(draw, draw > self.min_max[1])
            if len(draw) > 0:
                rtn = _np.append(rtn, draw)
                required = n - len(rtn)
        return _round_samples(rtn, round_to_decimals)


class Normal2D(PyNSNDistribution):

    def __init__(self, mu, sigma, correlation, max_radius=None):
        """Two dimensional normal distribution with optional radial cut-off

        Resulting multivariate distribution has the defined means, standard
        deviations and correlation (approx.) between the two variables

        Parameter:
        ----------
        mu:  tuple of to numbers (numeric, numeric)
        sigma: tuple of to numbers  (numeric, numeric)
        correlation: numeric
            the correlation between a and b
        max_radius: numerical (optional)
            cut-off criterion: maximal distance from the distribution center (mu)
        """
        super().__init__(min_max=(None, None))

        try:
            lmu = len(mu)
            lsigma = len(sigma)
        except:
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
        """variance covariance matrix"""
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
                r = np.hypot(draw[:,0], draw[:,1])
                draw = _np.delete(draw, r > self.max_radius, axis=0)

            if len(draw) > 0:
                # append
                if rtn is None:
                    rtn = draw
                else:
                    rtn = _np.append(rtn, draw, axis=0)
                required = n - len(rtn)

        return _round_samples(rtn, round_to_decimals)

    def as_dict(self):
        d = super().as_dict()
        d.update({"mu": self.mu.tolist(),
                  "sigma": self.sigma.tolist(),
                  "correlation": self.correlation,
                  "max_radius": self.max_radius})
        return d


class Beta(_PyNSNDistributionMuSigma):

    def __init__(self, mu=None, sigma=None, alpha=None, beta=None, min_max=None):
        """Beta distribution defined by the number range, min_max=(min, max),
         the mean and the standard deviation (std)

        Resulting distribution has the defined mean and std

        Parameter:
        ----------
        mu: numeric
        sigma: numeric
        min_max : tuple (numeric, numeric)
            the range of the distribution

        Note:
        -----
        Depending on the position of the mean in number range the
        distribution is left or right skewed.

        See Also:
        --------
        for calculated shape parameters [alpha, beta] see property
        `shape_parameter_beta`
        """
        if alpha is not None and beta is not None and (mu, sigma) == (None, None):
            mu, sigma = Beta._calc_mu_sigma(alpha, beta, min_max)
        elif mu is None or sigma is None or alpha is not None or beta is not None:
            raise TypeError("Either Mu & Sigma or Alpha & Beta have to specified.")
        super().__init__(mu=mu, sigma=sigma, min_max=min_max)

    def sample(self, n, round_to_decimals=None):
        if self.sigma is None or self.sigma == 0:
            return _np.array([self.mu] * n)

        alpha, beta = self.shape_parameter
        dist = _rng.generator.beta(a=alpha, b=beta, size=n)
        dist = (dist - _np.mean(dist)) / _np.std(dist)  # z values
        rtn = dist * self.sigma + self.mu
        return _round_samples(rtn, round_to_decimals)

    @property
    def shape_parameter(self):
        """Alpha (p) & beta (q) parameter for the beta distribution
        http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm

        Returns
        -------
        parameter: tuple
            shape parameter (alpha, beta) of the distribution

        """
        r = float(self.min_max[1] - self.min_max[0])
        m = (self.mu - self.min_max[0]) / r  # mean
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
    def _calc_mu_sigma(alpha, beta, min_max):
        a = float(alpha)
        b = float(beta)
        r = float((min_max[1] - min_max[0]))

        e = a / (a + b)
        mu = e * r + min_max[0]

        v = (a * b) / ((a + b) ** 2 * (a + b + 1))
        sigma = _np.sqrt(v) * r
        return mu, sigma
