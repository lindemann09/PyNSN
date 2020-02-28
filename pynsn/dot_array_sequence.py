"""
Dot Array Sequence
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from hashlib import md5 as _md5
import numpy as _np

from pynsn._lib import _misc as _misc
from pynsn._lib._dot_array import DotArray as _DotArray
from pynsn._lib._visual_features import Features as _Feat
from pynsn._lib import random_dot_array

class DASequence(object):

    def __init__(self):
        """ docu the use of numerosity_idx see get_array_numerosity
        dot array, please use append_dot_array and delete_array to modify the list

        """

        self.dot_arrays = []
        self.method = None
        self.error = None
        self.numerosity_idx = {}

    def append_dot_arrays(self, arr):
        if isinstance(arr, _DotArray):
            arr = [arr]
        self.dot_arrays.extend(arr)
        self.numerosity_idx = {da.features.numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def delete_dot_arrays(self, array_id):
        self.dot_arrays.pop(array_id)
        self.numerosity_idx = {da.features.numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def get_array(self, numerosity):
        """returns array with a particular numerosity"""

        try:
            return self.dot_arrays[self.numerosity_idx[numerosity]]
        except:
            return None

    @property
    def min_max_numerosity(self):
        return (self.dot_arrays[0].features.numerosity, self.dot_arrays[-1].features.numerosity)

    @property
    def hash(self):
        """meta hash of object ids"""

        m = _md5()
        for da in self.dot_arrays:
            m.update(da.hash.encode("UTF-8"))
        return m.hexdigest()

    def get_features_dict(self):
        """dictionary with arrays"""

        dicts = [x.features.get_features_dict() for x in self.dot_arrays]
        rtn = _misc.join_dict_list(dicts)
        rtn['sequence_id'] = [self.hash] * len(self.dot_arrays)  # all arrays have the same sequence ID
        return rtn

    def get_features_dataframe(self):
        from pandas import DataFrame, np
        d = self.get_features_dict()
        array = []
        for x in range(len(d["Hash"])):
            # transposing (numpy not possible, because of data types
            row = map(lambda k: d[k][x], d.keys())
            array.append(list(row))
        return DataFrame(array, columns=list(d.keys()))

    def get_numerosity_correlations(self):
        feat = self.get_features_dict()
        del feat['hash']
        feat_np = _np.round(_np.array(feat.values()).T, 2)
        cor = _np.corrcoef(feat_np, rowvar=False)
        cor = cor[0, :]
        names = feat.keys()
        rtn = {}
        for x in range(1, len(cor)):
            rtn[names[x]] = cor[x]
        return rtn

    def __str__(self):
        return self.get_csv()

    def get_csv(self, variable_names=True, colour_column=False,
                hash_column=True):

        rtn = ""
        tmp_var_names = variable_names

        for da in self.dot_arrays:
            rtn += da.get_csv(num_idx_column=True, hash_column=False,
                              variable_names=tmp_var_names,
                              colour_column=colour_column)
            tmp_var_names = False

        if hash_column:
            obj_id = self.hash
            rtn2 = ""
            tmp_var_names = variable_names
            for l in rtn.split("\n"):
                if tmp_var_names:
                    rtn2 += "hash," + l + "\n"
                    tmp_var_names = False
                elif len(l) > 0:
                    rtn2 += "{},{}\n".format(obj_id, l)
            return rtn2
        else:
            return rtn


def create(specs,
           match_feature,
           match_value,
           min_max_numerosity,
           round_decimals = None,
           source_number = None):  # todo could be an iterator
    """factory function

    Methods takes take , you might use make Process
        match_properties:
                continuous property or list of continuous properties to be match
                or None
     returns False is error occured (see self.error)
    """
    try:
        l = len(min_max_numerosity)
    except:
        l = 0
    if l != 2:
        raise ValueError("min_max_numerosity has to be a pair of (min, max)")

    min, max = sorted(min_max_numerosity)

    if source_number is None:
        if match_feature in [_Feat.SPARSITY]:
            source_number = min
        elif match_feature in [_Feat.FIELD_AREA, _Feat.COVERAGE] :
            source_number = max
        else:
            source_number = min + ((max - min)//2)

    _misc.check_feature_list(match_feature)

    if not isinstance(specs, random_dot_array.Specs):
        raise TypeError("Specs have to be random_dot_array.Specs, but not {}".format(
            type(specs).__name__))

    # keep field area
    if match_feature in list(_Feat.SPACE_FEATURES) + [_Feat.COVERAGE]:
        prefer_keeping_field_area = True
    else:
        prefer_keeping_field_area = False

    # make source stimulus
    if source_number is None:
        source_number = min + int((max - min)/2)
    source_da = random_dot_array.create(n_dots=source_number,
                                        specs=specs)
    source_da.match.match_feature(feature=match_feature, value=match_value)
    source_da.center_array()
    source_da.round(round_decimals)

    # matched deviants
    rtn = DASequence()
    rtn.method = match_feature

    # decreasing
    if min < source_number:
        tmp, error = _make_matched_deviants(
            reference_da=source_da,
            match_feature=match_feature,
            target_numerosity=min,
            round_decimals=round_decimals,
            prefer_keeping_field_area=prefer_keeping_field_area)

        rtn.append_dot_arrays(list(reversed(tmp)))
        if error is not None:
            rtn.error = error
    # source number
    rtn.append_dot_arrays(source_da)
    # increasing
    if max > source_number:
        tmp, error = _make_matched_deviants(
            reference_da=source_da,
            match_feature=match_feature,
            target_numerosity=max,
            round_decimals=round_decimals,
            prefer_keeping_field_area=prefer_keeping_field_area)
        rtn.append_dot_arrays(tmp)
        if error is not None:
            rtn.error = error

    return rtn

def _make_matched_deviants(reference_da, match_feature, target_numerosity,
                           round_decimals, prefer_keeping_field_area):
    """helper function. Do not use this method. Please use make"""



    if reference_da.features.numerosity == target_numerosity:
        change = 0
    elif reference_da.features.numerosity > target_numerosity:
        change = -1
    else:
        change = 1

    da = reference_da.copy()
    da_sequence = []

    error = None
    #print(match_props, target_numerosity)
    while True:
        try:
            da = da.number_deviant(change_numerosity=change,
                                   prefer_keeping_field_area=prefer_keeping_field_area)
        except:
            return [], "ERROR: Can't find the a make matched deviants"

        da.match.match_feature(feature=match_feature,
                                   reference_dot_array=reference_da)
        cnt = 0
        while True:
            cnt += 1
            ok, mesg = da.realign()
            if ok:
                break
            if cnt > 10:
                error = u"ERROR: realign, " + str(cnt) + ", " + str(da.features.numerosity)

        #print(da.feature.get_features_text())
        da.round(round_decimals)
        da_sequence.append(da)

        if error is not None or da.features.numerosity == target_numerosity:
            break

    return da_sequence, error

