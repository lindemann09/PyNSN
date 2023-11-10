import numpy as np
from numpy.typing import NDArray


def perimeter(sizes: NDArray) -> NDArray[np.float_]:
    """Ramanujan's second approximation of the ellipse perimeter"""
    s = np.atleast_2d(sizes)
    a = s[:, 0]
    b = s[:, 1]
    return np.pi * ((a+b) + (3*(a-b)**2) / (10*(a+b) + np.sqrt(a**2 + 14*a*b + b**2)))


def diameter(size: NDArray, theta: NDArray) -> NDArray[np.float_]:
    """Ellipse diameter at a certain angle

    Parameter
    ---------
        Size: NDarray
            2d array with (semi-majors, semi-minors)
        theta: float
            angle in radians
    """
    d = np.atleast_2d(size)
    return (d[:, 0] * d[:, 1]) / np.sqrt((d[:, 0] * np.sin(theta))**2
                                         + (d[:, 1] * np.cos(theta))**2)