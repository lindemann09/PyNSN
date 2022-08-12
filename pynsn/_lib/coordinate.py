from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from numpy.typing import NDArray, ArrayLike
import numpy as np


class Coordinate(object):
    __slots__ = ("_xy", )

    def __init__(self, xy: ArrayLike) -> None:
        self._xy = np.empty(2)
        self.xy = xy  # call setter

    def __repr__(self) -> str:
        return f"Coordinate(xy={self._xy})"

    def __add__(self, other: Coordinate) -> Coordinate:
        return Coordinate(self._xy + other._xy)

    def __sub__(self, other: Coordinate) -> Coordinate:
        return Coordinate(self._xy - other._xy)

    def __mul__(self, other: Coordinate) -> Coordinate:
        return Coordinate(np.multiply(self._xy, other._xy))

    def __div__(self, other: Coordinate) -> Coordinate:
        return Coordinate(np.divide(self._xy, other._xy))

    def __iadd__(self, other: Coordinate) -> Coordinate:
        self._xy = self._xy + other._xy
        return self

    def __isub__(self, other: Coordinate) -> Coordinate:
        self._xy = self.xy - other.xy  # type: ignore
        return self

    def __imul__(self, other: Coordinate) -> Coordinate:
        self._xy = np.multiply(self._xy, other._xy)
        return self

    def __idiv__(self, other: Coordinate) -> Coordinate:
        self._xy = np.divide(self._xy, other.xy)
        return self

    def __eq__(self, other: Coordinate) -> bool:
        return np.array_equal(self._xy, other._xy, equal_nan=True)

    def __ne__(self, other: Coordinate) -> bool:
        return not np.array_equal(self._xy, other._xy, equal_nan=True)

    @property
    def xy(self) -> NDArray[np.floating]:
        return self._xy  # type: ignore

    @xy.setter
    def xy(self, value: ArrayLike) -> None:
        value = np.asarray(value)
        if value.shape != (2,):
            raise ValueError("xy has be an iterable object with two elements")
        self._xy = value

    @property
    def x(self):
        return self._xy[0]

    @property
    def y(self):
        return self._xy[1]

    @property
    def rho(self) -> float:
        """rho (i.e., the radius) of the polar coordinate"""
        return np.hypot(self._xy[0], self._xy[1])

    @rho.setter
    def rho(self, value) -> None:
        self.polar = (value, self.theta)

    @property
    def theta(self) -> float:
        """theta (i.e., the angle) of the polar coordinate"""
        return np.arctan2(self._xy[1], self._xy[0])

    @theta.setter
    def theta(self, value) -> None:
        self.polar = (self.rho, value)

    @property
    def polar(self) -> NDArray:
        """Polar coordinate: rho (radius), theta (angle)"""
        return np.array([self.rho, self.theta])

    @polar.setter
    def polar(self, rad_ang) -> None:
        self._xy = np.array([rad_ang[0] * np.cos(rad_ang[1]),
                            rad_ang[0] * np.sin(rad_ang[1])])

    def distance(self, other: Coordinate) -> float:
        """Euclidean distance to the another Coordinate."""
        d_xy = self._xy - other.xy
        return np.hypot(d_xy[0], d_xy[1])
