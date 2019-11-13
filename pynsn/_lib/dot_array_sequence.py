"""
Dot Array Sequence
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from hashlib import md5 as _md5
import numpy as _np

from copy import copy as _copy
from . import misc as _misc
from . import visual_features as _vf
from ._dot_array import DotArray as _DotArray

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
        self.numerosity_idx = {da.feature.numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def delete_dot_arrays(self, array_id):
        self.dot_arrays.pop(array_id)
        self.numerosity_idx = {da.feature.numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def get_array(self, numerosity):
        """returns array with a particular numerosity"""

        try:
            return self.dot_arrays[self.numerosity_idx[numerosity]]
        except:
            return None

    @property
    def min_max_numerosity(self):
        return (self.dot_arrays[0].feature.numerosity, self.dot_arrays[-1].feature.numerosity)

    @property
    def hash(self):
        """meta hash of object ids"""

        m = _md5()
        for da in self.dot_arrays:
            m.update(da.hash.encode("UTF-8"))
        return m.hexdigest()

    def get_features_dict(self):  # todo search for get_features!
        """dictionary with arrays"""

        dicts = [x.feature.get_features_dict() for x in self.dot_arrays]
        rtn = _misc.join_dict_list(dicts)
        rtn['hash'] = [self.hash] * len(self.dot_arrays)  # all arrays have the
        # same ID
        return rtn

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
                picture_column=False, hash_column=True):

        rtn = ""
        tmp_var_names = variable_names

        for da in self.dot_arrays:
            rtn += da.get_csv(num_idx_column=True, hash_column=False,
                              variable_names=tmp_var_names,
                              colour_column=colour_column,
                              picture_column=picture_column)
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


def create(reference_dot_array,
           match_properties,
           min_max_numerosity,
           extra_space,  #  later fitting of convex hull and density might
           # result in enlarged arrays TODO rename variable
           center_array=True,
           logger=None):  # todo could be an iterator
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

    if match_properties is None:
        match_properties = []
    elif not isinstance(match_properties, (tuple, list)):
        match_properties = [match_properties]
    _vf.check_feature_list(match_properties, check_set_value=False)

    if not isinstance(reference_dot_array, _DotArray):
        raise TypeError("Reference_dot_array has to be DotArray, but not {}".format(
            type(reference_dot_array).__name__))


    # copy and change values to match this stimulus
    match_props = []
    prefer_keeping_field_area = False
    for m in match_properties:
        m = _copy(m)
        m.adapt_value(reference_dot_array)
        match_props.append(m)
        if isinstance(m, _vf.LogSpacing().dependencies) or \
                (isinstance(m, _vf.Coverage) and
                 m.match_ratio_fieldarea2totalarea < 1):
            prefer_keeping_field_area = True
            break

    # adjust reference (basically centering)
    reference_da = reference_dot_array.copy()
    reference_da.target_array_radius += (extra_space // 2)  # add extra space
    if center_array:
        reference_da._xy -= reference_da.center_of_outer_positions
        reference_da.set_array_modified()

    # matched deviants
    rtn = DASequence()
    rtn.method = match_props

    min, max = sorted(min_max_numerosity)
    # decreasing
    if min < reference_dot_array.feature.numerosity:
        da_sequence, error = _make_matched_deviants(
            reference_da=reference_da,
            match_props=match_props,
            target_numerosity=min,
            prefer_keeping_field_area=prefer_keeping_field_area)
        rtn.append_dot_arrays(list(reversed(da_sequence)))
        if error is not None:
            rtn.error = error
    # reference
    rtn.append_dot_arrays(reference_da)
    # increasing
    if max > reference_dot_array.feature.numerosity:
        da_sequence, error = _make_matched_deviants(
            reference_da=reference_da,
            match_props=match_props,
            target_numerosity=max,
            prefer_keeping_field_area=prefer_keeping_field_area)
        rtn.append_dot_arrays(da_sequence)
        if error is not None:
            rtn.error = error

    if logger is not None:
        from ._logging import LogFile # to avoid circular import
        if not isinstance(logger, LogFile):
            raise TypeError("logger has to be None or a GeneratorLogger, and not {}".format(
                type(logger).__name__))

        logger.log(rtn)

    return rtn

def _make_matched_deviants(reference_da, match_props, target_numerosity,
                           prefer_keeping_field_area):  # TODO center array OK?
    """helper function. Do not use this method. Please use make"""

    if reference_da.feature.numerosity == target_numerosity:
        change = 0
    elif reference_da.feature.numerosity > target_numerosity:
        change = -1
    else:
        change = 1

    da = reference_da.copy()
    da_sequence = []

    error = None
    while True:
        try:
            da = da.number_deviant(change_numerosity=change,
                               prefer_keeping_field_area=prefer_keeping_field_area)
        except:
            return [], "ERROR: Can't find the a make matched deviants"
        _vf.check_feature_list(match_props)
        for feat in match_props:
            da.match(feat)

        cnt = 0
        while True:
            cnt += 1
            ok, mesg = da.realign()
            if ok:
                break
            if cnt > 10:
                error = u"ERROR: realign, " + str(cnt) + ", " + str(da.feature.numerosity)

        da_sequence.append(da)

        if error is not None or da.feature.numerosity == target_numerosity:
            break

    return da_sequence, error

