import sys as _sys
import typing as _tp

import numpy as _np
import numpy.typing as _ntp
from ._collections import CollectionStimulusPairs
from .. import rnd as _rnd
from .._stimulus.properties import VisProp as _VisProp
from .._stimulus.properties import ensure_vis_prop as _ensure_vis_prop
from .. import fit as _stim_fit


def property_ratio_correlation(collection: CollectionStimulusPairs,
                               distr: _tp.Union[_rnd.AbstractUnivarDistr, _rnd.Abstract2dDistr],
                               prop_a: _tp.Union[str, _VisProp],
                               prop_b: _tp.Union[None, str, _VisProp] = None,
                               max_corr: float = 0.01,
                               feedback: bool = True) -> _tp.Union[_tp.Tuple[float, float], float]:

    prop_a = _ensure_vis_prop(prop_a)
    prop_b = _ensure_vis_prop(prop_b)  # FIXME None not yet supported

    num_ratios = collection.property_ratios(_VisProp.N).to_numpy()
    num_ratios = collection.property_dataframe()["N"].to_numpy()
    rnd_values, target_correlations = _get_rnd_target_values(
        num_ratios, distr=distr, prop_a=prop_a, prop_b=prop_b,
        max_corr=max_corr)

    # print(rnd_values)
    # print(collection.property_ratios([prop_a, prop_b]))
    n = len(collection.pairs)
    for i, sp in enumerate(collection.pairs):
        before = sp.property_difference
        if feedback:
            _sys.stdout.write(
                f"fitting {i+1}/{n} {sp.name}                 \r")
        _stim_fit.property_ratio(sp, prop_a, rnd_values[i, 0])
        if isinstance(prop_b, _VisProp):
            _stim_fit.property_ratio(sp, prop_b, rnd_values[i, 1])
    if feedback:
        print(" "*70)

    # print(collection.property_ratios([prop_a, prop_b]))
    return target_correlations


def property_difference_correlation(collection: CollectionStimulusPairs,
                                    distr: _tp.Union[_rnd.AbstractUnivarDistr, _rnd.Abstract2dDistr],
                                    prop_a: _tp.Union[str, _VisProp],
                                    prop_b: _tp.Union[None,
                                                      str, _VisProp] = None,
                                    max_corr: float = 0.01,
                                    feedback: bool = True) -> _tp.Union[_tp.Tuple[float, float], float]:

    prop_a = _ensure_vis_prop(prop_a)
    prop_b = _ensure_vis_prop(prop_b)  # FIXME None not yet supported

    num_dist = collection.property_differences(_VisProp.N).to_numpy()
    rnd_values, target_correlations = _get_rnd_target_values(
        num_dist, distr=distr, prop_a=prop_a, prop_b=prop_b,
        max_corr=max_corr)

    n = len(collection.pairs)
    for i, sp in enumerate(collection.pairs):
        if feedback:
            _sys.stdout.write(
                f"fitting {i+1}/{n} {sp.name}                 \r")
        _stim_fit.property_difference(sp, prop_a, rnd_values[i, 0])
        if isinstance(prop_b, _VisProp):
            _stim_fit.property_difference(sp, prop_b, rnd_values[i, 1])
    if feedback:
        print(" "*70)

    return target_correlations

# helper


def _get_rnd_target_values(number_list: _ntp.NDArray,
                           distr: _tp.Union[_rnd.AbstractUnivarDistr, _rnd.Abstract2dDistr],
                           prop_a: _VisProp,
                           prop_b: _tp.Optional[_VisProp] = None,
                           max_corr=0.01):
    if isinstance(prop_b, _VisProp):
        if prop_a.is_dependent_from(prop_b):
            raise ValueError(f"'{prop_a.name}' and '{prop_b.name}' depend" +
                             " on each other and can't be varied independently")
        if isinstance(distr, _rnd.Abstract2dDistr):
            return _modify_2d_distributions(distr,
                                            number_list=number_list,
                                            max_corr=max_corr)
        else:
            raise ValueError("distr has to be a 2 dimensional distribution")
    else:
        if isinstance(distr, _rnd.AbstractUnivarDistr):
            return _modify_distributions(distr,
                                         number_list=number_list,
                                         max_corr=max_corr)
        else:
            raise ValueError("distr has to be a univariate distribution")


def _modify_2d_distributions(dist: _rnd.Abstract2dDistr,
                             number_list: _ntp.NDArray,
                             max_corr=0.01) -> _tp.Tuple[_ntp.NDArray[_np.float64], _tp.Tuple[float, float]]:

    n = len(number_list)
    while True:
        values = dist.sample(n)
        r1 = _np.corrcoef(number_list, values[:, 0])[0, 1]
        if _np.abs(r1) <= max_corr:
            r2 = _np.corrcoef(number_list, values[:, 1])[0, 1]
            if _np.abs(r2) <= max_corr:
                return values, (float(r1), float(r2))


def _modify_distributions(dist: _rnd.AbstractUnivarDistr,
                          number_list: _ntp.NDArray,
                          max_corr=0.01,) -> _tp.Tuple[_ntp.NDArray[_np.float64], float]:

    n = len(number_list)
    while True:
        values = dist.sample(n)
        r = _np.corrcoef(number_list, values)[0, 1]
        if _np.abs(r) <= max_corr:
            return _np.atleast_2d(values).T, float(r)
