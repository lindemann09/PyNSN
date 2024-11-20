from __future__ import annotations

import json
from pathlib import Path
from typing import Union

import numpy as np

from .. import defaults
from .nsn_stimulus import NSNStimulus
from .properties import VisProp, ensure_vis_prop


class NSNStimulusPair():

    def __init__(self,
                 stim_a: NSNStimulus,
                 stim_b: NSNStimulus,
                 name: str = "no_name") -> None:
        self.stim_a = stim_a
        self.stim_b = stim_b
        self.name = name

    def to_json(self,
                path: Union[None, str, Path] = None,
                indent: int = 2, tabular: bool = True):
        """Save the StimulusPair as folder with json files"""

        a = self.stim_a.to_json(indent=indent, tabular=tabular)
        b = self.stim_b.to_json(indent=indent, tabular=tabular)
        c = '{"name":"' + self.name + '"}'
        json_str = f"[{c},\n{a},\n{b}]"
        if isinstance(path, (Path, str)):
            with open(path, "w", encoding=defaults.FILE_ENCODING) as fl:
                fl.write(json_str)
        return json_str

    @staticmethod
    def from_json(path: Union[str, Path]) -> NSNStimulusPair:
        """Load StimulusPair from json file

        see `to_json`
        """
        path = Path(path)
        if not path.is_file():
            raise RuntimeError(f"Can't find {path}.")

        with open(path, 'r', encoding=defaults.FILE_ENCODING) as fl:
            c, a, b = json.load(fl)

        return NSNStimulusPair(stim_a=NSNStimulus.from_dict(a),
                               stim_b=NSNStimulus.from_dict(b),
                               name=c["name"])

    def property_difference(self, prop: Union[str, VisProp]) -> np.float64:
        """difference of property `prop` between stimulus A & B"""
        prop = ensure_vis_prop(prop)
        rtn = self.stim_a.properties.get(
            prop) - self.stim_b.properties.get(prop)
        return np.float64(rtn)

    def property_ratio(self, prop: Union[str, VisProp]) -> np.float64:
        """ratio of property `prop` between stimulus A & B"""
        prop = ensure_vis_prop(prop)
        rtn = self.stim_a.properties.get(
            prop) / self.stim_b.properties.get(prop)
        return np.float64(rtn)
