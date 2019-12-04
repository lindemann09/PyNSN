import os
import json
import gzip
from ._lib._dot_array import DotArray
from .dot_array_sequence import DASequence

def load(json_file_name, zipped=False):

    if not os.path.exists(json_file_name):
        new_flname = json_file_name + ".json"
        if os.path.exists(new_flname):
            json_file_name = new_flname
        else:
            raise RuntimeError("File '{}' does not exists.".format(json_file_name))

    if zipped:
        fl = gzip.open(json_file_name, 'r')
    else:
        fl = open(json_file_name, 'r')

    rtn = DotArrayArchive()
    rtn.dict = json.load(fl)
    fl.close()
    return rtn

class DotArrayArchive(object):

    def __init__(self):
        self.dict = {}

    def add(self, dot_array):

        if isinstance(dot_array, DotArray):
            self.dict[dot_array.hash] = dot_array.as_dict()
        elif isinstance(dot_array, DASequence):
            hash_list = list(map(lambda x: x.hash, dot_array.dot_arrays))
            self.dict[dot_array.hash] = {"sequence" : hash_list}
            for da in dot_array.dot_arrays:
                self.add(da)
        else:
            RuntimeError("dot_array has to be a pynsn.DotArray or a "
                         "pynsn.dot_array_sequence.DASequence")

    @property
    def array_ids(self):
        rtn = []
        for k,v in self.dict.items():
            if "xy" in v:
                rtn.append(k)
        return rtn

    @property
    def sequence_ids(self):
        rtn = []
        for k,v in self.dict.items():
            if "sequence" in v:
                rtn.append(k)
        return rtn

    def get_dot_array(self, id):
        try:
            d = self.dict[id]
        except:
            return None
        if "xy" not in d:
            return None # It is a DASequene

        rtn = DotArray(0,0)
        rtn.read_from_dict(d)
        return rtn

    def get_da_sequence(self, id):
        try:
            d = self.dict[id]
        except:
            return None
        if "sequence" not in d:
            return None

        tmp = []
        for id in d["sequence"]:
            tmp.append(self.get_dot_array(id))
            if tmp[-1] is None:
                raise RuntimeError("Can't find dot_array {}".format(id))
        rtn = DASequence()
        rtn.append_dot_arrays(tmp)
        return rtn

    def _all_features(self):
        """helper function returns array with all features and varnames"""
        array = []
        feat = {}
        for id in self.array_ids:
            feat = self.get_dot_array(id).features.get_features_dict()
            array.append(list(feat.values()))

        varnames = map(lambda x:x.replace(" ", "_"), feat.keys())
        return array, list(varnames)

    def features_csv(self, delimiter =","):

        array, varnames = self._all_features()
        rtn = delimiter.join(varnames)
        for row in array:
            rtn += "\n" + delimiter.join(map(lambda x:str(x), row))
        return rtn

    def features_dataframe(self):
        """returns pandas dataframe with all features"""
        from pandas import DataFrame
        array, varnames = self._all_features()
        return DataFrame(array, columns=varnames)

    def save(self, json_file_name, indent=None, zipped=False):
        if not (json_file_name.endwith(".json") or
                json_file_name.endwith(".json.gz")):
            json_file_name += ".json"
            if zipped:
                json_file_name += ".gz"

        if zipped:
            fl = gzip.open(json_file_name, 'w')
        else:
            fl = open(json_file_name, 'w')

        json.dump(self.dict, fl, indent=indent)
        fl.close()
