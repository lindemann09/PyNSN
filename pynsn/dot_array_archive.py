import json
from ._lib._dot_array import DotArray
from .dot_array_sequence import DASequence

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

    def get(self):
        return self.dict

    def save(self, filename, indent=None):
        with open(filename, 'w') as fl:
            json.dump(self.dict, fl, indent=indent)