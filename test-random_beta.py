from matplotlib.pyplot import *
import random as _random
import numpy as np
from nosynum import *

def random_norm(a, b, sigma = None):
    """normal distributed random value between a and b """
    if sigma is None:
        sigma = (b-a)/6.0
    r = _random.normalvariate(mu = a + (b-a)/2.0, sigma = sigma)
    if r<a or r>b:
        return random_norm(a=a, b=b)
    return r

def coin_flip(head_bias = 0.5):
    """Return randomly True (head) or False (tail).

    Parameters
    ----------
    head_bias : numeric, optional
        bias in favor of head (default=0.5, fair coin)

    Returns
    -------
    rnd : bool

    """

    if head_bias<0 or head_bias>1:
        raise RuntimeError("Head bias must be between 0 and 1!")

    return _random.random()<= head_bias




number_range=(180, 200)
small_mean_relative_to_interval = 0.25 # percent in interval
large_mean_relative_to_interval = 0.75
std_in_percent_of_range = 0.15

r = number_range[1] - number_range[0]
mean_small = number_range[0] + r * small_mean_relative_to_interval
mean_large = number_range[0] + r * large_mean_relative_to_interval
std = r * std_in_percent_of_range


#parameter = shape_parameter_beta(number_range, mean_small, std)
#print "parameter", parameter
#x = [random_beta(number_range, parameter) for _ in range(10000)]
#hist(x, bins=100)
#print (np.mean(x), np.std(x))

#parameter = shape_parameter_beta(number_range, mean_large, std)
parameter = [2, 2]

#x = [random_beta(number_range, parameter) for _ in range(10000)]
x = [int(coin_flip(0.7)) for _ in range(100000)]

hist(x, bins=100)
print (np.mean(x), np.std(x))


show()
