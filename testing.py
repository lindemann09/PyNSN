import pynsn
from pynsn import random
import matplotlib.pyplot as plt

d = random.Normal2D(mu=(0, 0), sigma=(10, 10),
                    xy_minmax=(-10, None, 10, 20),
                    correlation=0.5)
d.pyplot_samples()

plt.show()
