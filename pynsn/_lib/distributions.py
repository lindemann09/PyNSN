import random as _random
# NOTE: do not use numpy.random, because it produces identical numbers for
# different threads
import numpy as _np

_random.seed()



class _PyNSNDistribution(object):

    def __init__(self, min_max):
        if not isinstance(min_max, (list, tuple)) or len(min_max) != 2 or \
                (None not in min_max and min_max[0] > min_max[1]):
            raise TypeError("min_max {} ".format(min_max) + \
                            "has to be a tuple of two values (a, b) with a <= b.")
        self.min_max = min_max

    def as_dict(self):
        return {"type": type(self).__name__,
                "min_max": self.min_max}

    def _cutoff_outside_range(self, np_vector):
        # helper function that cuts off values outside min_max range
        if self.min_max[0] is not None:
            np_vector = _np.delete(np_vector, np_vector < self.min_max[0])
        if self.min_max[1] is not None:
            np_vector = _np.delete(np_vector, np_vector > self.min_max[1])
        return np_vector

    def sample(self, n):
        pass

class _PyNSNDistributionMuSigma(_PyNSNDistribution):

    def __init__(self, mu, sigma, min_max):
        super().__init__(min_max)
        self.mu = mu
        self.sigma = abs(sigma)
        if (min_max[0] is not None and mu <= min_max[0]) or \
                (min_max[1] is not None and mu >= min_max[1]):
            txt = "mean ({}) has to be inside the defined min_max range ({})".format(
                mu, min_max)
            raise RuntimeError(txt)

    def as_dict(self):
        d = super().as_dict()
        d.update({"mu" : self.mu,
                  "sigma": self.sigma})
        return d


class Laplace(_PyNSNDistribution):
    """
    """
    def __init__(self, min_max):
        """Laplace distribution defined by the number range, min_max=(min, max)

        Parameter:
        ----------
        min_max : tuple (numeric, numeric)
            the range of the distribution
        """
        super().__init__(min_max)

    def sample(self, n):
        dist = _np.array([_random.random() for _ in range(n)])
        r = float(self.min_max[1] - self.min_max[0])
        return self.min_max[0] + dist * r


class Normal(_PyNSNDistributionMuSigma):

    def __init__(self,  mu, sigma, min_max=None):
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

    def sample(self, n):
        rtn = _np.array([])
        required = n
        while required>0:
            draw = _np.array([_random.normalvariate(self.mu, self.sigma) \
                          for _ in range(required)])
            rtn = self._cutoff_outside_range(_np.append(rtn, draw))
            required = n - len(rtn)
        return rtn


class Beta(_PyNSNDistributionMuSigma):

    def __init__(self, mu, sigma, min_max):
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
        super().__init__(mu=mu, sigma=sigma, min_max=min_max)

    def sample(self, n):
        if self.sigma is None or self.sigma == 0:
            return _np.array([self.mu] * n)

        alpha, beta = self.shape_parameter
        dist = _np.array([_random.betavariate(alpha=alpha, beta=beta) \
                         for _ in range(n)])
        dist = (dist - _np.mean(dist)) / _np.std(dist) # z values
        return dist * self.sigma + self.mu

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
        mean = (self.mu - self.min_max[0]) / r
        std = self.sigma / r
        x = (mean * (1 - mean) / std ** 2 - 1)
        return x * mean, (1 - mean) * x

