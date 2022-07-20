"""
Dot Array Sequence
"""
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

#FIXME Depricated model currently BROKEN
# FIXME BROKEN MODULE

from hashlib import md5 as _md5
import numpy as _np

from .._lib import misc as _misc
from .._lib import DotArray
from ..visual_properties import flags, fit


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
        self.dotextend(arr)
        self.numerosity_idx = {da._properties.numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def delete_dot_arrays(self, array_id):
        self.dotpop(array_id)
        self.numerosity_idx = {da._properties.numerosity: idx for idx, da in enumerate(self.dot_arrays)}

    def get_array(self, numerosity):
        """returns array with a particular numerosity"""

        try:
            return self.dot_arrays[self.numerosity_idx[numerosity]]
        except:
            return None

    @property
    def min_max_numerosity(self):
        return self.dot_arrays[0].properties.numerosity, \
               self.dot_arrays[-1].properties.numerosity

    @property
    def hash(self):
        """meta hash of object ids"""

        m = _md5()
        for da in self.dot_arrays:
            m.update(da.hash.encode("UTF-8"))
        return m.hexdigest()

    def get_properties_dict(self):
        """dictionary with arrays

        Examples
        --------
        making a pandas dataframe with aa properties
        >>> d = my_da_sequence.as_dict()
        >>> array = []
        >>> for x in range(len(d["Hash"])):
        >>>    row = map(lambda k: d[k][x], d.keys())
        >>>    array.append(list(row))
        >>> return pandas.DataFrame(array, columns=list(d.keys()))
        """
        dicts = [x._properties.as_dict() for x in self.dot_arrays]
        rtn = _misc.join_dict_list(dicts)
        rtn['sequence_id'] = [self.hash] * len(self.dot_arrays)  # all arrays have the same _sequence ID
        return rtn


    def get_numerosity_correlations(self):
        feat = self.get_properties_dict()
        del feat['hash']
        prop_np = _np.round(_np.asarray(feat.values()).T, 2)
        cor = _np.corrcoef(prop_np, rowvar=False)
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
            rtn += da.csv(num_idx_column=True, hash_column=False,
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
           adapt_property,
           adapt_value,
           min_max_numerosity,
           round_decimals = None,
           source_number = None):  # todo could be an iterator
    """factory function

    Methods takes take , you might use make Process
        adapt_properties:
                continuous property or list of continuous properties to be adapt
                or None
     returns False is error occurred (see self.error)
    """
    try:
        l = len(min_max_numerosity)
    except:
        l = 0
    if l != 2:
        raise ValueError("min_max_numerosity has to be a pair of (min, max)")

    min_, max_ = sorted(min_max_numerosity)

    if source_number is None:
        if adapt_property in [flags.SPARSITY]:
            source_number = min_
        elif adapt_property in [flags.FIELD_AREA, flags.COVERAGE] :
            source_number = max_
        else:
            source_number = min_ + ((max_ - min_)//2)

    check_property_list(adapt_property)

    if not isinstance(specs, DotArraySpecs):
        raise TypeError("Specs has to be of type DotArraySpecs, but not {}".format(
            type(specs).__name__))

    # keep field area
    if adapt_property in list(flags.SPACE_FEATURES) + [flags.COVERAGE]:
        prefer_keeping_field_area = True
    else:
        prefer_keeping_field_area = False

    # make source image
    if source_number is None:
        source_number = min_ + int((max_ - min_)/2)
    source_da = random_array.create(n_objects=source_number,
                                    specs=specs)
    source_da = fit.visual_property(source_da, property_flag=adapt_property, value=adapt_value)
    source_da.center_array()
    source_da.mod_round_values(round_decimals)

    # adapted deviants
    rtn = DASequence()
    rtn.method = adapt_property

    # decreasing
    if min_ < source_number:
        tmp, error = _make_adapted_deviants(
            reference_da=source_da,
            adapt_property=adapt_property,
            target_numerosity=min_,
            round_decimals=round_decimals,
            prefer_preserve_field_area=prefer_keeping_field_area)

        rtn.append_dot_arrays(list(reversed(tmp)))
        if error is not None:
            rtn.error = error
    # source number
    rtn.append_dot_arrays(source_da)
    # increasing
    if max_ > source_number:
        tmp, error = _make_adapted_deviants(
            reference_da=source_da,
            adapt_property=adapt_property,
            target_numerosity=max_,
            round_decimals=round_decimals,
            prefer_preserve_field_area=prefer_keeping_field_area)
        rtn.append_dot_arrays(tmp)
        if error is not None:
            rtn.error = error

    return rtn

def _make_adapted_deviants(reference_da, adapt_property, target_numerosity,
                           round_decimals, prefer_preserve_field_area):
    """helper function. Do not use this method. Please use make"""



    if reference_da._properties.numerosity == target_numerosity:
        change = 0
    elif reference_da._properties.numerosity > target_numerosity:
        change = -1
    else:
        change = 1

    da = reference_da.copy()
    da_sequence = []

    error = None
    #print(adapt_props, target_numerosity)
    while True:
        try:
            da = da.get_number_deviant(change_numerosity=change,
                                       preserve_field_area=prefer_preserve_field_area)
        except:
            return [], "ERROR: Can't find the a make adapted deviants"

        da = fit.visual_property(da, property_flag=adapt_property,
                                 value=reference_da.properties.get(adapt_property))
        cnt = 0
        while True:
            cnt += 1
            ok, mesg = da.realign()
            if ok:
                break
            if cnt > 10:
                error = u"ERROR: realign, " + str(cnt) + ", " + str(da._properties.numerosity)

        da.mod_round_values(round_decimals)
        da_sequence.append(da)

        if error is not None or da._properties.numerosity == target_numerosity:
            break

    return da_sequence, error


def check_property_list(feature_list):
    """helper function
    raises TypeError or Runtime errors if checks fail
    * type check
    * dependency check
    """

    size_occured = ""
    space_occured = ""
    error = "Incompatible properties to adapt: {} & {}"

    if not isinstance(feature_list, (tuple, list)):
        feature_list = [feature_list]

    for x in feature_list:
<<<<<<< HEAD:pynsn/_sequence/dot_array_sequence.py
        if x not in _Feat.ALL_FEATURES:
            raise ValueError("Parameter is not a continuous feature or a " + \
=======
        if x not in flags:
            raise ValueError("Parameter is not a continuous propertyor a " + \
>>>>>>> devel:pynsn/tools/_dot_array_sequence.py
                            "list of continuous properties")
            # continious property or visual prop

        if x in flags.SIZE_FEATURES:
            if len(size_occured)>0:
                raise ValueError(error.format(x, size_occured))
            else:
                size_occured = x

        if x in flags.SPACE_FEATURES:
            if len(space_occured)>0:
                raise RuntimeError(error.format(x, space_occured))
            else:
                space_occured = x

