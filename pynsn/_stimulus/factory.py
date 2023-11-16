"""
"""
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import deepcopy
from typing import Optional, Union

from .. import defaults
from .._misc import formated_json
from .._shapes import Dot, Ellipse, PolygonShape, Rectangle
from .._shapes.abc_shapes import AbstractShape
from ..random._random_shape import AbstractRndShape
from .nsn_stimulus import NSNStimulus
from .stimulus_colours import StimulusColours

# todo incremental_random_dot_array


class StimulusFactory(object):
    """Factory class for creating Non-Symbolic Number Stimulus
    """

    def __init__(self,
                 target_area_shape: Union[Dot, Rectangle,
                                          Ellipse, PolygonShape, None] = None,
                 min_distance: int = defaults.MIN_DISTANCE,
                 min_distance_target_area: int = defaults.MIN_DISTANCE,
                 stimulus_colours: Optional[StimulusColours] = None,
                 ignore_overlaps: bool = False,
                 max_iterations: Optional[int] = None
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
        self._stim = NSNStimulus(target_area_shape=target_area_shape,
                                 min_distance=min_distance,
                                 min_distance_target_area=min_distance_target_area)
        if stimulus_colours is not None:
            self._stim.colours = stimulus_colours
            if target_area_shape.colour.value is not None:
                self._stim.colours.target_area = target_area_shape.colour

        if max_iterations is None:
            self.max_iterations = defaults.MAX_ITERATIONS
        else:
            self.max_iterations = max_iterations

        self.ignore_overlaps = ignore_overlaps
        self._shapes = []

    @property
    def base_stimulus(self) -> NSNStimulus:
        return self._stim

    @base_stimulus.setter
    def base_stimulus(self, stim: NSNStimulus) -> None:
        assert isinstance(stim, NSNStimulus)
        self._stim = stim

    def clear(self):
        self._shapes = []

    def add(self, shape: Union[AbstractShape, AbstractRndShape], n: int = 1):
        assert (isinstance(shape, AbstractShape) or
                isinstance(shape, AbstractRndShape))
        self._shapes.append((n, shape))

    def get(self) -> NSNStimulus:
        rtn = deepcopy(self._stim)
        for n, s in self._shapes:
            rtn.add_shapes(ref_object=s, n=n, random_position=True,
                           ignore_overlaps=self.ignore_overlaps, max_iterations=self.max_iterations)
        return rtn

    def todict(self) -> dict:
        """dict representation of the stimulus factory"""
        rtn = {}
        rtn.update({"type": self.__class__.__name__})
        base_dict = self._stim.todict()
        del base_dict["shape_table"]
        del base_dict["hash"]
        rtn.update(base_dict)
        for x, s in enumerate(self._shapes):
            d = {"n": s[0]}
            d.update(s[1].todict())
            rtn.update({f"shape{x}": d})
        return rtn

    def tojson(self, filename: str = "", indent: int = 2) -> str:
        """json representation of the stimulus factory"""
        json_str = formated_json(self.todict(), indent=indent)
        if filename:
            with open(filename, "w", encoding="utf-8") as fl:
                fl.write(json_str)
        return json_str
