"""
Draw a random number from a beta dirstribution
"""
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import random

random.seed()

def random_beta(number_range, parameter):
	"""Draw from beta distribution defined by the
	number_range [a,b] and the shape parameters [alpha, beta]

	Parameter:
	----------
	number_range : tuple (numeric, numeric)
		the range of the distribution
	parameter: tuple
		shape parameter (alpha, beta) of the distribution
		see parameter function()

	Note:
	-----
	Depending on the position of the mean in number range the
	distribution is left or right skewed.

	"""

	return number_range[0] + (number_range[1] - number_range[0])\
			* random.betavariate(alpha=parameter[0], beta=parameter[1])

def shape_parameter_beta(number_range, mean, std):
	"""Returns alpha (p) & beta (q) parameter for the beta distribution
	http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm

	Parameter
	---------
	number_range : tuple (numeric, numeric)
		the range of the distribution
	mean : numeric
		the distribution mean
	std : numeric
		 the distribution standard deviation

	Returns
	-------
	parameter: tuple
		shape parameter (alpha, beta) of the distribution

	"""

	if mean<=number_range[0] or mean>=number_range[1] or \
		number_range[0]>=number_range[1]:
		raise RuntimeError("Mean has to be inside the defined number range")
	f = float(number_range[1] - number_range[0])
	mean = (mean-number_range[0])/f
	std = (std)/f
	x = ( mean*(1-mean) / std**2 - 1 )
	return (x*mean, (1-mean)*x)
