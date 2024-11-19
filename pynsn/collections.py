from __future__ import annotations

import typing as _tp
from pathlib import Path as _Path

import pandas as _pd

from . import NSNStimulus, NSNStimulusPair, _misc

ListNSNStimuli = _tp.List[NSNStimulus]
ListNSNStimPairs = _tp.List[NSNStimulusPair]


class CollectionStimulusPairs():
    """set of NSNNumPairs"""

    def __init__(self, lst: _tp.Union[None, ListNSNStimPairs] = None) -> None:

        if isinstance(lst, _tp.List):
            for x in lst:
                if not isinstance(x, NSNStimulusPair):
                    raise RuntimeError(
                        f"lst must be a list of NSNStimulusPairs and not {type(x)}")
            self.pairs = lst
        else:
            self.pairs: ListNSNStimPairs = []
        self._prop_a = _pd.DataFrame()
        self._prop_b = _pd.DataFrame()

    def append(self, stim_a: NSNStimulus, stim_b: NSNStimulus,
               name: str = "no_name"):
        """append two stimulus to the collection
        """

        self.pairs.append(NSNStimulusPair(stim_a, stim_b, name))
        self._prop_a = _pd.DataFrame()
        self._prop_b = _pd.DataFrame()

    def to_json(self, folder: _tp.Union[str, _Path]):
        """Save the collection as json files organized in subfolder"""
        for x in self.pairs:
            x.to_json(folder)

    @staticmethod
    def from_json(folder: _tp.Union[str, _Path]) -> CollectionStimulusPairs:
        """Load collection from subfolders with json files

        see `to_json`
        """
        folder = _Path(folder)
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
            pa["stim"] = "a"
            pa["name"] = x.name

            pb = x.stim_a.properties.to_dict(True)
            pb["stim"] = "b"
            pb["name"] = x.name
            a.append(pa)
            b.append(pb)

        self._prop_a = _pd.DataFrame(_misc.dict_of_arrays(a))
        self._prop_b = _pd.DataFrame(_misc.dict_of_arrays(b))

    def property_differences(self,
                             columns: _tp.Union[str, _tp.List[str]]) -> _pd.DataFrame:
        """differences of property `prop` between stimuli A & B"""

        if len(self._prop_a) != len(self.pairs):
            self.calc_properties()
        return self._prop_a[columns] - self._prop_b[columns]

    def property_ratios(self,
                        columns: _tp.Union[str, _tp.List[str]]) -> _pd.DataFrame:
        """ratios of property `prop` between stimuli A & B"""

        if len(self._prop_a) != len(self.pairs):
            self.calc_properties()
        return self._prop_a[columns] / self._prop_b[columns]
