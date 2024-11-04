
import timeit
import numpy as np
import matplotlib.pyplot as plt
from pynsn import rnd
import pynsn
a = range(10)
index = (8, 4, 2)
new = []
for c, s in enumerate(a):
    if c not in index:
        new.append(s)
print(a)
print(new)
exit()


start = timeit.timeit()

d = rnd.Uniform2D(x_minmax=(-100, 100),
                     y_minmax=(20, 50))

s = d.sample(1)
print(s)

# d.pyplot_samples()

# plt.show()

# d = random.Uniform2D(x_minmax=(50, 55), y_minmax=(-100, 100),
#                     radial_radius=6)
