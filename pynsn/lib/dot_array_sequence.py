"""
Dot Array Sequence
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from hashlib import md5
import numpy as np

from copy import copy
from . import misc, visual_features
from .dot_array import DotArray

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
        if isinstance(arr, DotArray):
            arr = [arr]
        self.dot_arrays.extend(arr)
        self.numerosity_idx = {da.feature_numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def delete_dot_arrays(self, array_id):
        self.dot_arrays.pop(array_id)
        self.numerosity_idx = {da.feature_numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def get_array(self, numerosity):
        """returns array with a particular numerosity"""

        try:
            return self.dot_arrays[self.numerosity_idx[numerosity]]
        except:
            return None

    @property
    def min_max_numerosity(self):
        return (self.dot_arrays[0].feature_numerosity, self.dot_arrays[-1].feature_numerosity)

    @property
    def object_id(self):
        """meta hash of object ids"""

        m = md5()
        for da in self.dot_arrays:
            m.update(da.object_id.encode("UTF-8"))
        return m.hexdigest()[:DotArray.OBJECT_ID_LENGTH]

    def get_features_dict(self):  # todo search for get_features!
        """dictionary with arrays"""

        dicts = [x.get_features_dict() for x in self.dot_arrays]
        rtn = misc.join_dict_list(dicts)
        rtn['object_id'] = [self.object_id] * len(self.dot_arrays)  # all arrays have the same ID
        return rtn

    def get_numerosity_correlations(self):
        feat = self.get_features_dict()
        del feat['object_id']
        feat_np = np.round(np.array(feat.values()).T, 2)
        cor = np.corrcoef(feat_np, rowvar=False)
        cor = cor[0, :]
        names = feat.keys()
        rtn = {}
        for x in range(1, len(cor)):
            rtn[names[x]] = cor[x]
        return rtn

    def __str__(self):
        return self.get_csv()

    def get_csv(self, num_format="%7.2f", variable_names=True, colour_column=False,
                picture_column=False, object_id_column=True):

        rtn = ""
        tmp_var_names = variable_names
        for da in self.dot_arrays:
            rtn += da.get_csv(num_idx_column=True, object_id_column=False,
                              variable_names=tmp_var_names,
                              num_format=num_format, colour_column=colour_column,
                              picture_column=picture_column)
            tmp_var_names = False

        if object_id_column:
            obj_id = self.object_id
            rtn2 = ""
            tmp_var_names = variable_names
            for l in rtn.split("\n"):
                if tmp_var_names:
                    rtn2 += "object_id," + l + "\n"
                    tmp_var_names = False
                elif len(l) > 0:
                    rtn2 += "{},{}\n".format(obj_id, l)
            return rtn2
        else:
            return rtn


def generate_da_sequence(reference_dot_array,
                  match_properties,
                  min_max_numerosity,
                  extra_space,  # fitting convex hull and density might result in enlarged arrays
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
    visual_features.check_feature_list(match_properties, check_set_value=False)

    if not isinstance(reference_dot_array, DotArray):
        raise TypeError("Reference_dot_array has to be DotArray, but not {}".format(
            type(reference_dot_array).__name__))


    # copy and change values to match this stimulus
    match_props = []
    prefer_keeping_field_area = False
    for m in match_properties:
        m = copy(m)
        m.adapt_value(reference_dot_array)
        match_props.append(m)
        if isinstance(m, visual_features.LogSpacing().dependencies) or \
                (isinstance(m, visual_features.Coverage) and m.match_ratio_fieldarea2totalarea < 1):
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
    if min < reference_dot_array.feature_numerosity:
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
    if max > reference_dot_array.feature_numerosity:
        da_sequence, error = _make_matched_deviants(
            reference_da=reference_da,
            match_props=match_props,
            target_numerosity=max,
            prefer_keeping_field_area=prefer_keeping_field_area)
        rtn.append_dot_arrays(da_sequence)
        if error is not None:
            rtn.error = error

    if logger is not None:
        from .logging import LogFile # to avoid circular import
        if not isinstance(logger, LogFile):
            raise TypeError("logger has to be None or a GeneratorLogger, and not {}".format(
                type(logger).__name__))

        logger.log(rtn)

    return rtn

def _make_matched_deviants(reference_da, match_props, target_numerosity,
                           prefer_keeping_field_area):  # TODO center array OK?
    """helper function. Do not use this method. Please use make"""

    if reference_da.feature_numerosity == target_numerosity:
        change = 0
    elif reference_da.feature_numerosity > target_numerosity:
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
        if len(match_props) > 0:
            da.match(match_props, realign=False, center_array=False)

        cnt = 0
        while True:
            cnt += 1
            ok, mesg = da.realign()
            if ok:
                break
            if cnt > 10:
                error = u"ERROR: realign, " + str(cnt) + ", " + str(da.feature_numerosity)

        da_sequence.append(da)

        if error is not None or da.feature_numerosity == target_numerosity:
            break

    return da_sequence, error

