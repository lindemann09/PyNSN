import pynsn
from pynsn import random
import matplotlib.pyplot as plt
import numpy as np
import timeit
start = timeit.timeit()

d = random.Uniform2D(x_minmax=(-100, 100),
                     y_minmax=(20, 50))

s = d.sample(100000)
print(np.min(s, axis=0))
print(np.max(s, axis=0))

# d.pyplot_samples()

# plt.show()

# d = random.Uniform2D(x_minmax=(50, 55), y_minmax=(-100, 100),
#                     radial_radius=6)
