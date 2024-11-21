from __future__ import annotations

import gzip
import json
import typing as tp
from pathlib import Path

import pandas as pd

from .. import _misc, defaults
from .._stimulus import NSNStimulus
from .._stimulus.properties import VPList, ensure_vp
from ._abc_coll import AbstractCollection, ListNSNStimuli

# FIXME not yet tested
# TODO NSNStimuli have no names yet


class CollectionStimuli(AbstractCollection):
    """Collection of NSNNumPairs"""

    def __init__(self, lst: None | ListNSNStimuli = None) -> None:

        if isinstance(lst, tp.List):
            for x in lst:  # type check
                if not isinstance(x, NSNStimulus):
                    raise RuntimeError(
                        f"lst must be a list of NSNStimulus and not {type(x)}")
            self.stimuli = lst
        else:
            self.stimuli: ListNSNStimuli = []
        self._prop_df = pd.DataFrame()

    def append(self, stim: NSNStimulus):
        """append stimulus to the collection
        """

        self.stimuli.append(stim)
        self.reset_properties()

    def save(self, path: tp.Union[str, Path], zipped: bool = True):
        """Save the collection as json files organized in subfolder"""

        json_str = "["
        for x in self.stimuli:
            s = x.to_json(path=None, indent=2, tabular=True)
            json_str += f"\n{s},"
        json_str = json_str[:-1] + "\n]"

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if zipped:
            with gzip.open(path, "wt", encoding=defaults.FILE_ENCODING) as fl:
                fl.write(json_str)
        else:
            with open(path, "w", encoding=defaults.FILE_ENCODING) as fl:
                fl.write(json_str)

    @staticmethod
    def load(path: tp.Union[str, Path], zipped: bool = True) -> CollectionStimuli:
        """Load collection from folder with json files

        see `to_json`
        """
        path = Path(path)
        if not path.is_file():
            raise RuntimeError(f"Can't find {path}.")
        if zipped:
            with gzip.open(path, 'rt', encoding=defaults.FILE_ENCODING) as fl:
                dicts = json.load(fl)
        else:
            with open(path, 'r', encoding=defaults.FILE_ENCODING) as fl:
                dicts = json.load(fl)

        rtn = CollectionStimuli()
        for d in dicts:
            rtn.append(NSNStimulus.from_dict(d))
        return rtn

    def reset_properties(self):
        """resets visual properties

        If the array `CollectionStimuli.stimuli` have been changed directly,
        the method needs to be called to ensure get valid property data
        """
        self._prop_df = pd.DataFrame()

    def _calc_properties(self):
        """re-calculate all`visual properties,
        should be called by private methods after `reset_properties`"""

        a = []
        for x in self.stimuli:
            pa = x.properties.to_dict(True)
            a.append(pa)

        self._prop_df = pd.DataFrame(_misc.dict_of_arrays(a))

    def property_dataframe(self) -> pd.DataFrame:
        """Dataframe with all properties of the two stimuli (`a` & `b`)

        The dataframe contains additional the name of the pairs
        """

        if len(self._prop_df) != len(self.stimuli):
            self._calc_properties()

        rtn = self._prop_df.copy()
        rtn["names"] = [x.name for x in self.stimuli]
        return rtn

    def corr(self, properties: None | VPList = None) -> pd.DataFrame:
        """Pairwise Pearson correlation between  of the visual properties"""
        if len(self._prop_df) != len(self.stimuli):
            self._calc_properties()

        df = self._prop_df
        if properties is not None:
            prop_names = [ensure_vp(p).name for p in properties]
            df = df[prop_names]
        return df.corr()
