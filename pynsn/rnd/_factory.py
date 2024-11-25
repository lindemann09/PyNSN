"""
"""

from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .. import defaults
from .._misc import formated_json
from .._shapes import Dot, Ellipse, PolygonShape, Rectangle
from .._shapes.abc_shapes import AbstractShape
from .._stimulus.nsn_stimulus import NSNStimulus
from .._stimulus.stimulus_colours import StimulusColours
from ._random_shape import AbstractRndShape

# TODO  incremental_random_dot_array

FactoryShapesList = List[Tuple[int, AbstractShape | AbstractRndShape]]


class StimulusFactory(object):
    """Factory class for creating Non-Symbolic Number Stimulus
    """

    def __init__(self,
                 target_area_shape: Dot | Rectangle | Ellipse | PolygonShape | None = None,
                 min_distance: int = defaults.MIN_DISTANCE,
                 min_distance_target_area: int = defaults.MIN_DISTANCE,
                 stimulus_colours: StimulusColours | None = None,
                 ignore_overlaps: bool = False,
                 max_iterations: int | None = None
                 ) -> None:
        """

        Args:
            target_shape_or_stim: _description_
            min_distance: _description_. Defaults to defaults.MIN_DISTANCE.
            min_distance_target_area: _description_. Defaults to defaults.MIN_DISTANCE.
            stimulus_colours: _description_. Defaults to None.
            ignore_overlaps: _description_. Defaults to False.
            max_iterations: _description_. Defaults to None.

        """
        if not isinstance(target_area_shape, AbstractShape):
            # default target area
            target_area_shape = Rectangle(size=(400, 400))
        self._base = NSNStimulus(target_area_shape=target_area_shape,
                                 min_distance=min_distance,
                                 min_distance_target_area=min_distance_target_area)
        if stimulus_colours is not None:
            self._base.colours = stimulus_colours
            if target_area_shape.colour.value is not None:
                self._base.colours.target_area = target_area_shape.colour

        if max_iterations is None:
            self.max_iterations = defaults.MAX_ITERATIONS
        else:
            self.max_iterations = max_iterations

        self.ignore_overlaps = ignore_overlaps
        self._shapes: FactoryShapesList = []

    @property
    def base_stimulus(self) -> NSNStimulus:
        return self._base

    @base_stimulus.setter
    def base_stimulus(self, stim: NSNStimulus) -> None:
        assert isinstance(stim, NSNStimulus)
        self._base = stim

    def shapes(self) -> FactoryShapesList:
        return self._shapes

    def shapes_clear(self):
        self._shapes: FactoryShapesList = []

    def shape_add(self, shape: AbstractShape | AbstractRndShape, n: int = 1):
        assert (isinstance(shape, AbstractShape) or
                isinstance(shape, AbstractRndShape))
        self._shapes.append((n, shape))

    def create(self, name: str | None = None,
               extra_distance_target_area: int | None = None,
               ) -> NSNStimulus:
        """Create a random stimulus

        Parameters
        ----------
        name : str | None, optional
            _description_, by default None
        extra_distance_target_area : int | None, optional
            extra minimum distance to the boarder of the target area, by default None.

            The random process considers this plus the stimulus property `min_distance_target_area`
            as defined while instantiation of StimulusFactory


        Returns
        -------
        NSNStimulus
        """

        rtn = deepcopy(self._base)
        if extra_distance_target_area is not None:
            rtn.target_area.min_dist_boarder += extra_distance_target_area

        if name is not None:
            rtn.name = name
        for n, s in self._shapes:
            rtn.shape_add_random_pos(ref_object=s, n=n,
                                     ignore_overlaps=self.ignore_overlaps, max_iterations=self.max_iterations)
        return rtn

    def to_dict(self) -> dict:
        """dict representation of the stimulus factory"""
        rtn = {}
        rtn.update({"type": self.__class__.__name__})
        base_dict = self._base.to_dict(tabular=False)

        for key in ("shape_array", "hash", "name"):
            base_dict.pop(key, None)  # avoid KeyError

        rtn.update(base_dict)
        for x, s in enumerate(self._shapes):
            d = {"n": s[0]}
            d.update(s[1].to_dict())
            rtn.update({f"shape{x}": d})
        return rtn

    def to_json(self,
                filename: str | Path | None = None,
                indent: int = 2,
                decimals: None | int = None) -> str:
        """json representation of the stimulus factory"""
        json_str = formated_json(
            self.to_dict(), indent=indent, decimals=decimals)
        if isinstance(filename, (Path, str)):
            with open(filename, "w", encoding=defaults.FILE_ENCODING) as fl:
                fl.write(json_str)
        return json_str

    @staticmethod
    def from_dict(the_dict: Dict[str, Any]) -> StimulusFactory:
        raise NotImplementedError()  # FIXME
