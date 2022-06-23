import random
# NOTE: do not use numpy.random, because it produces identical numbers for
# different threads
import numpy as np
from .misc import numpy_round2

random.seed()

class PyNSNDistribution(object):

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
            np_vector = np.delete(np_vector, np_vector < self.min_max[0])
        if self.min_max[1] is not None:
            np_vector = np.delete(np_vector, np_vector > self.min_max[1])
        return np_vector

    def sample(self, n, round_to_decimals=False):
        return np.array([0] * n)

    def pyplot_samples(self, n=100000):

        try:
            from matplotlib.pyplot import hist
        except:
            raise ImportError("To use pyplot, please install matplotlib.")

        return hist(self.sample(n=n), bins=100)[2]


class _PyNSNDistributionMuSigma(PyNSNDistribution):

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


class Uniform(PyNSNDistribution):
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

    def sample(self, n, round_to_decimals=None):
        dist = np.array([random.random() for _ in range(n)])
        r = float(self.min_max[1] - self.min_max[0])
        rtn = self.min_max[0] + dist * r
        if round_to_decimals is not None:
            return numpy_round2(rtn, decimals=round_to_decimals)
        else:
            return rtn


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

    def sample(self, n, round_to_decimals=None):
        rtn = np.array([])
        required = n
        while required>0:
            draw = np.array([random.normalvariate(self.mu, self.sigma) \
                             for _ in range(required)])
            rtn = self._cutoff_outside_range(np.append(rtn, draw))
            required = n - len(rtn)

        if round_to_decimals is not None:
            return numpy_round2(rtn, decimals=round_to_decimals)
        else:
            return rtn


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
            return np.array([self.mu] * n)

        alpha, beta = self.shape_parameter
        dist = np.array([random.betavariate(alpha=alpha, beta=beta) \
                         for _ in range(n)])
        dist = (dist - np.mean(dist)) / np.std(dist) # z values
        rtn = dist * self.sigma + self.mu

        if round_to_decimals is not None:
            return numpy_round2(rtn, decimals=round_to_decimals)
        else:
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
        r = float(self.min_max[1] - self.min_max[0])
        m = (self.mu - self.min_max[0]) / r # mean
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
        mu =  e * r + min_max[0]

        v = (a*b) / ((a+b)**2 * (a+b+1))
        sigma = np.sqrt(v) * r
        return mu, sigma

