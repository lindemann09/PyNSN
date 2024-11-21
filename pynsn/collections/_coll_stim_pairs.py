from __future__ import annotations

import gzip
import json
import typing as tp
from pathlib import Path

import pandas as pd

from .. import _misc, defaults
from .._stimulus import NSNStimulus, NSNStimulusPair
from .._stimulus.properties import ensure_vp, VPOrList
from ._abc_coll import AbstractCollection, ListNSNStimPairs
from ._coll_stim import CollectionStimuli


class CollectionStimulusPairs(AbstractCollection):
    """Collection of NSNNumPairs"""

    def __init__(self, lst: tp.Union[None, ListNSNStimPairs] = None) -> None:

        if isinstance(lst, tp.List):
            for x in lst:
                if not isinstance(x, NSNStimulusPair):
                    raise RuntimeError(
                        f"lst must be a list of NSNStimulusPairs and not {type(x)}")
            self.pairs = lst
        else:
            self.pairs: ListNSNStimPairs = []
        self._prop_df_a = pd.DataFrame()
        self._prop_df_b = pd.DataFrame()

    def append(self, stim_a: NSNStimulus, stim_b: NSNStimulus,
               name: str = "no_name"):
        """append two stimulus to the collection
        """

        self.pairs.append(NSNStimulusPair(stim_a, stim_b, name))
        self.reset_properties_dataframe()

    def save(self, path: tp.Union[str, Path], zipped: bool = True):
        """Save the collection as json files organized in subfolder"""

        json_str = "["
        for x in self.pairs:
            s = x.to_json(path=None, indent=2, tabular=True)
            json_str += "\n" + s[1:-1] + ","
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
    def load(path: tp.Union[str, Path], zipped: bool = True) -> CollectionStimulusPairs:
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

        rtn = CollectionStimulusPairs()
        while len(dicts) > 0:
            c = dicts.pop(0)
            a = dicts.pop(0)
            b = dicts.pop(0)
            rtn.append(NSNStimulus.from_dict(a),
                       NSNStimulus.from_dict(b),
                       name=c["name"])
        return rtn

    def reset_properties_dataframe(self):
        """reset dataframe of visual properties

        If the array `CollectionStimulusPairs.pairs` have been changed directly,
        the method needs to be called to ensure get valid property data
        """
        self._prop_df_a = pd.DataFrame()
        self._prop_df_b = pd.DataFrame()

    def _calc_properties(self):
        """re-calculate all`visual properties,
        should be called by private methods after `reset_properties`"""

        a = []
        b = []
        for x in self.pairs:
            pa = x.stim_a.properties.to_dict(True)
            pb = x.stim_b.properties.to_dict(True)
            a.append(pa)
            b.append(pb)

        self._prop_df_a = pd.DataFrame(_misc.dict_of_arrays(a))
        self._prop_df_b = pd.DataFrame(_misc.dict_of_arrays(b))

    def property_dataframe(self) -> pd.DataFrame:
        """Dataframe with all properties of the two stimuli (`a` & `b`)

        The dataframe contains additional the name of the pairs
        """

        if len(self._prop_df_a) != len(self.pairs):
            self._calc_properties()

        names = [x.name for x in self.pairs]
        a = self._prop_df_a.copy()
        a["stim"] = "a"
        a["name"] = names
        b = self._prop_df_b.copy()
        b["stim"] = "b"
        b["name"] = names
        return pd.concat([a, b])

    def property_differences(self, props:  None | VPOrList = None) -> pd.DataFrame:
        """differences of properties between stimuli A & B

        optionally specify `props`
        """

        if len(self._prop_df_a) != len(self.pairs):
            self._calc_properties()

        if props is None:
            return self._prop_df_a - self._prop_df_b
        elif isinstance(props, tp.List):
            cols = [ensure_vp(p).name for p in props]
        else:
            cols = ensure_vp(props).name

        return self._prop_df_a[cols] - self._prop_df_b[cols]

    def property_ratios(self, props:  None | VPOrList = None) -> pd.DataFrame:
        """ratios of properties between stimuli A & B

        optionally specify `props`
        """

        if len(self._prop_df_a) != len(self.pairs):
            self._calc_properties()

        if props is None:
            return self._prop_df_a / self._prop_df_b
        elif isinstance(props, tp.List):
            cols = [ensure_vp(p).name for p in props]
        else:
            cols = ensure_vp(props).name

        return self._prop_df_a[cols] / self._prop_df_b[cols]

    def collection_stimuli(self, stim_id: str) -> CollectionStimuli:
        """
        stim_id : {'a', 'b'}
        """

        rtn = CollectionStimuli()
        for pair in self.pairs:
            if stim_id == "a":
                rtn.append(pair.stim_a)
            elif stim_id == "b":
                rtn.append(pair.stim_b)
            else:
                raise ValueError(f"Unknown stimulus id '{stim_id}'. " +
                                 "Must be either 'a' or 'b'")
        return rtn
