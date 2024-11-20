from __future__ import annotations

import typing as tp
from pathlib import Path

import pandas as pd

from .. import _misc
from .._stimulus import NSNStimulus, NSNStimulusPair
from .._stimulus.properties import (VisProp, VisPropList, ensure_vis_prop,
                                    ensure_vis_prop_list)

ListNSNStimuli = tp.List[NSNStimulus]
ListNSNStimPairs = tp.List[NSNStimulusPair]


OptionalVisPropList = None | VisProp | str | VisPropList


class CollectionStimulusPairs():
    """set of NSNNumPairs"""

    def __init__(self, lst: tp.Union[None, ListNSNStimPairs] = None) -> None:

        if isinstance(lst, tp.List):
            for x in lst:
                if not isinstance(x, NSNStimulusPair):
                    raise RuntimeError(
                        f"lst must be a list of NSNStimulusPairs and not {type(x)}")
            self.pairs = lst
        else:
            self.pairs: ListNSNStimPairs = []
        self._prop_a = pd.DataFrame()
        self._prop_b = pd.DataFrame()

    def append(self, stim_a: NSNStimulus, stim_b: NSNStimulus,
               name: str = "no_name"):
        """append two stimulus to the collection
        """

        self.pairs.append(NSNStimulusPair(stim_a, stim_b, name))
        self._prop_a = pd.DataFrame()
        self._prop_b = pd.DataFrame()

    def to_json(self, folder: tp.Union[str, Path]):
        """Save the collection as json files organized in subfolder"""
        for x in self.pairs:
            x.to_json(folder)

    @staticmethod
    def from_json(folder: tp.Union[str, Path]) -> CollectionStimulusPairs:
        """Load collection from subfolders with json files

        see `to_json`
        """
        folder = Path(folder)
        if not folder.is_dir():
            raise RuntimeError(
                f"Can't load from {folder}. It's not a directory")

        rtn = CollectionStimulusPairs()
        for d in folder.iterdir():
            if d.is_dir():
                rtn.pairs.append(NSNStimulusPair.from_json(d))

        return rtn

    def calc_properties(self):
        """calculate all visual properties

        If the array `CollectionStimulusPairs.pairs` have been changed directly,
        the method needs to be called to ensure get valid property data
        """

        a = []
        b = []
        for x in self.pairs:
            pa = x.stim_b.properties.to_dict(True)
            pb = x.stim_a.properties.to_dict(True)
            a.append(pa)
            b.append(pb)

        self._prop_a = pd.DataFrame(_misc.dict_of_arrays(a))
        self._prop_b = pd.DataFrame(_misc.dict_of_arrays(b))

    def property_dataframe(self) -> pd.DataFrame:

        if len(self._prop_a) != len(self.pairs):
            self.calc_properties()

        names = [x.name for x in self.pairs]
        a = self._prop_a.copy()
        a["stim"] = "a"
        a["name"] = names
        b = self._prop_b.copy()
        b["stim"] = "b"
        b["name"] = names
        return pd.concat([a, b])

    def property_differences(self, props: OptionalVisPropList = None) -> pd.DataFrame:
        """differences of properties between stimuli A & B

        optionally specify `props`
        """

        if len(self._prop_a) != len(self.pairs):
            self.calc_properties()

        if props is None:
            return self._prop_a - self._prop_b
        elif isinstance(props, tp.List):
            cols = [p.as_string() for p in ensure_vis_prop_list(props)]
        else:
            cols = ensure_vis_prop(props).as_string()

        return self._prop_a[cols] - self._prop_b[cols]

    def property_ratios(self, props: OptionalVisPropList = None) -> pd.DataFrame:
        """ratios of properties  between stimuli A & B

        optionally specify `props`
        """

        if len(self._prop_a) != len(self.pairs):
            self.calc_properties()

        if props is None:
            return self._prop_a / self._prop_b
        elif isinstance(props, tp.List):
            cols = [p.as_string() for p in ensure_vis_prop_list(props)]
        else:
            cols = ensure_vis_prop(props).as_string()

        return self._prop_a[cols] / self._prop_b[cols]
