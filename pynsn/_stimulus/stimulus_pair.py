from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
from .nsn_stimulus import NSNStimulus
from .properties import VisProp


class NSNStimulusPair():

    def __init__(self,
                 stim_a: NSNStimulus,
                 stim_b: NSNStimulus,
                 name: str = "no_name") -> None:
        self.stim_a = stim_a
        self.stim_b = stim_b
        self.name = name

    def to_json(self, folder: Union[str, Path],
                indent: int = 2, tabular: bool = True):
        """Save the StimulusPair as folder with json files"""

        p = Path(folder, self.name)
        p.mkdir(parents=True, exist_ok=True)
        self.stim_a.to_json(Path(p, "a.json"), indent, tabular=tabular)
        self.stim_b.to_json(Path(p, "b.json"), indent, tabular=tabular)

    @staticmethod
    def from_json(path: Union[str, Path]) -> NSNStimulusPair:
        """Load StimulusPair from json files

        see `to_json`
        """
        path = Path(path)
        if not path.is_dir():
            raise RuntimeError(f"Can't load from {path}. It's not a directory")

        return NSNStimulusPair(stim_a=NSNStimulus.from_json(Path(path, "a.json")),
                               stim_b=NSNStimulus.from_json(
                                   Path(path, "b.json")),
                               name=path.stem)

    def property_difference(self, prop: VisProp) -> np.float64:
        """difference of property `prop` between stimulus A & B"""

        rtn = self.stim_a.properties.get(
            prop) - self.stim_b.properties.get(prop)
        return np.float64(rtn)

    def property_ratio(self, prop: VisProp) -> np.float64:
        """ratio of property `prop` between stimulus A & B"""

        rtn = self.stim_a.properties.get(
            prop) / self.stim_b.properties.get(prop)
        return np.float64(rtn)
