#FIXME Depricated model currently BROKEN

import os
import json
import gzip

<<<<<<< HEAD:pynsn/_sequence/dot_array_archive.py
from .._nsn.dot_array import DotArray
from pynsn._sequence.dot_array_sequence import DASequence
=======
from .._lib import DotArray
from ._dot_array_sequence import DASequence
>>>>>>> devel:pynsn/tools/_dot_array_archive.py

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

    rtn = DotArraySequenceArchive()
    rtn.dict = json.load(fl)
    fl.close()
    return rtn

class DotArraySequenceArchive(object):
    ## should actually use new sqlite achrive
    
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
            RuntimeError("nsn has to be a pynsn.DotArray or a "
                         "pynsn.dot_array_sequence.DASequence")

    def remove(self, id):
        try:
            self.dict.pop(id)
        except:
            return False
        return True

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
                raise RuntimeError("Can't find nsn {}".format(id))
        rtn = DASequence()
        rtn.append_dot_arrays(tmp)
        return rtn

    def all_properties(self):
        """eturns array with all properties and varnames

        Examples
        --------
        could be used to make a pandas dataframe

        >>> array, varnames = my_dot_array_archive.all_properties()
        >>> pandas.DataFrame(array, columns=varnames)
        """
        array = []
        feat = {}
        for id in self.array_ids:
            feat = self.get_dot_array(id)._properties.as_dict()
            array.append(list(feat.values()))

        varnames = map(lambda x:x.replace(" ", "_"), feat.keys())
        return array, list(varnames)

    def properties_csv(self, delimiter =","):

        array, varnames = self.all_properties()
        rtn = delimiter.join(varnames)
        for row in array:
            rtn += "\n" + delimiter.join(map(lambda x:str(x), row))
        return rtn

    def save(self, json_file_name, indent=None, zipped=False):

        if zipped:
            if not json_file_name.endswith(".json.gz"):
                if json_file_name.endswith(".json"):
                    json_file_name += ".gz"
                else:
                    json_file_name += ".json.gz"
            fl = gzip.open(json_file_name, "wb")

        else:
            if not json_file_name.endswith(".json"):
                json_file_name += ".json"
            fl = open(json_file_name, 'wb')

        fl.write(json.dumps(self.dict, indent=indent).encode("utf-8"))
        fl.close()

