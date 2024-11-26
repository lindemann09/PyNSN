import sys as _sys
from typing import Tuple as _Tuple

import numpy as _np
import numpy.typing as _ntp
from scipy.stats import linregress as _linregress
import pandas as _pd

from .. import fit as _stim_fit
from ..rnd._distributions import AbstractUnivarDistr as _AbstractUnivarDistr
from ..rnd._distributions_2d import Abstract2dDistr as _Abstract2dDistr
from .._stimulus.properties import VP as _VP
from .._stimulus.properties import ensure_vp as _ensure_vp
from ._coll_stim_pairs import CollectionStimuli as _CollectionStimuli
from ._coll_stim_pairs import CollectionStimulusPairs as _CollectionStimulusPairs


def property_correlation(stimuli: _CollectionStimuli,
                         distr: _AbstractUnivarDistr,
                         prop_a: str | _VP,
                         prop_b: None | str | _VP = None,
                         max_corr: float = 0.01,
                         feedback: bool = True) -> _Tuple[float, float] | float:
    prop_a = _ensure_vp(prop_a)
    if prop_b is not None:
        prop_b = _ensure_vp(prop_b)

    nums = _np.array([s.properties.numerosity for s in stimuli.stimuli])
    rnd_values, target_correlations = _check_prop_and_rnd_target_values(
        nums, distr=distr, prop_a=prop_a, prop_b=prop_b, max_corr=max_corr)

    n = len(stimuli.stimuli)
    for i, sp in enumerate(stimuli.stimuli):
        if feedback:
            _sys.stdout.write(
                f"fitting {i+1}/{n} {sp.name}                 \r")
        _stim_fit.property_adapt(sp, prop_a, rnd_values[i, 0])

        if isinstance(prop_b, _VP):
            _stim_fit.property_adapt(sp, prop_b, rnd_values[i, 1])

    if feedback:
        print(" "*70)

    stimuli.reset_properties()
    return target_correlations


def _minimize_correlation(ref_var:_pd.Series,
                          tar_var: _pd.Series,
                          omega: float = 1) -> _Tuple[_pd.Series, float]:
    n = len(tar_var)
    if len(ref_var) != n:
        raise ValueError("ref_var and tar_var have to have the same length!")

    cdata = _np.array([ref_var, tar_var]).T
    means = _np.mean(cdata, axis=0)
    for i in range(cdata.shape[1]):
        cdata[:, i] = cdata[:, i] - means[i]
    reg = _linregress(cdata[:, 0], cdata[:, 1])

    dat = tar_var - omega * reg.slope * cdata[:, 0]  # type: ignore
    reg = _linregress(cdata[:, 0], dat)
    return dat, float(reg.rvalue)  # type: ignore


def minimize_correlation_property_ratio(pairs: _CollectionStimulusPairs,
                                        prop_a: str | _VP,
                                        prop_b: None | str | _VP = None,
                                        omega: float = 1, # FIXME not needed
                                        adapt_stim: str = "both",
                                        feedback: bool = True) -> _Tuple[float, float] | float:
    prop_a = _ensure_vp(prop_a)
    props = pairs.property_ratios([_VP.N, prop_a])
    pa, target_correlations = _minimize_correlation(
        ref_var=props.loc[:, "N"], tar_var=props.loc[:, prop_a.name], omega=omega)

    if prop_b is not None:
        prop_b = _ensure_vp(prop_b)
        tmp = pairs.property_ratios([prop_b])
        pb, rb = _minimize_correlation(ref_var=props.loc[:, "N"],
                                       tar_var=tmp.loc[:, prop_b.name],
                                       omega=omega)
        target_correlations = (target_correlations, rb)
    else:
        pb = None

    # fitting
    n = len(pairs.pairs)
    for i, sp in enumerate(pairs.pairs):
        if feedback:
            _sys.stdout.write(
                f"fitting {i+1}/{n} {sp.name}                 \r")
        _stim_fit.property_ratio(sp, prop_a, pa[i], adapt_stim=adapt_stim)

        if isinstance(pb, _pd.Series):
            _stim_fit.property_ratio(
                sp, prop_b, pb[i], adapt_stim=adapt_stim)  # type: ignore

    if feedback:
        print(" "*70)

    pairs.reset_properties()
    return target_correlations


def property_ratio_correlation(pairs: _CollectionStimulusPairs,
                               distr: _AbstractUnivarDistr | _Abstract2dDistr,
                               prop_a: str | _VP,
                               prop_b: None | str | _VP = None,
                               max_corr: float = 0.01,
                               adapt_stim: str = "both",
                               feedback: bool = True) -> _Tuple[float, float] | float:

    prop_a = _ensure_vp(prop_a)
    if prop_b is not None:
        prop_b = _ensure_vp(prop_b)

    num_ratios = pairs.property_ratios(_VP.N).to_numpy()
    rnd_values, target_correlations = _check_prop_and_rnd_target_values(
        num_ratios, distr=distr, prop_a=prop_a, prop_b=prop_b, max_corr=max_corr)

    n = len(pairs.pairs)
    for i, sp in enumerate(pairs.pairs):
        if feedback:
            _sys.stdout.write(
                f"fitting {i+1}/{n} {sp.name}                 \r")
        _stim_fit.property_ratio(
            sp, prop_a, rnd_values[i, 0], adapt_stim=adapt_stim)

        if isinstance(prop_b, _VP):
            _stim_fit.property_ratio(
                sp, prop_b, rnd_values[i, 1], adapt_stim=adapt_stim)

    if feedback:
        print(" "*70)

    pairs.reset_properties()
    return target_correlations


def property_difference_correlation(pairs: _CollectionStimulusPairs,
                                    distr: _AbstractUnivarDistr | _Abstract2dDistr,
                                    prop_a: str | _VP,
                                    prop_b: None | str | _VP = None,
                                    max_corr: float = 0.01,
                                    feedback: bool = True) -> _Tuple[float, float] | float:

    prop_a = _ensure_vp(prop_a)
    if prop_b is not None:
        prop_b = _ensure_vp(prop_b)

    num_dist = pairs.property_differences(_VP.N).to_numpy()
    rnd_values, target_correlations = _check_prop_and_rnd_target_values(
        num_dist, distr=distr, prop_a=prop_a, prop_b=prop_b,
        max_corr=max_corr)

    n = len(pairs.pairs)
    for i, sp in enumerate(pairs.pairs):
        if feedback:
            _sys.stdout.write(
                f"fitting {i+1}/{n} {sp.name}                 \r")
        _stim_fit.property_difference(sp, prop_a, rnd_values[i, 0])
        if isinstance(prop_b, _VP):
            _stim_fit.property_difference(sp, prop_b, rnd_values[i, 1])
    if feedback:
        print(" "*70)

    pairs.reset_properties()

    return target_correlations

# helper


def _check_prop_and_rnd_target_values(number_list: _ntp.NDArray,
                                      distr: _AbstractUnivarDistr | _Abstract2dDistr,
                                      prop_a: _VP,
                                      prop_b: None | _VP = None,
                                      max_corr=0.01):
    if isinstance(prop_b, _VP):
        if prop_a.is_dependent_from(prop_b):
            raise ValueError(f"'{prop_a.name}' and '{prop_b.name}' depend" +
                             " on each other and can't be varied independently")
        if isinstance(distr, _Abstract2dDistr):
            return __modify_2d_distributions(distr,
                                             number_list=number_list,
                                             max_corr=max_corr)
        else:
            raise ValueError("distr has to be a 2 dimensional distribution," +
                             " if two properties should be fitted")
    else:
        if isinstance(distr, _AbstractUnivarDistr):
            return __modify_distributions(distr,
                                          number_list=number_list,
                                          max_corr=max_corr)
        else:
            raise ValueError("distr has to be a univariate distribution," +
                             " if one property should be fitted")


def __modify_2d_distributions(distr: _Abstract2dDistr,
                              number_list: _ntp.NDArray,
                              max_corr=0.01) -> _Tuple[_ntp.NDArray[_np.float64], _Tuple[float, float]]:

    n = len(number_list)
    while True:
        values = distr.sample(n)
        r1 = _np.corrcoef(number_list, values[:, 0])[0, 1]
        if _np.abs(r1) <= max_corr:
            r2 = _np.corrcoef(number_list, values[:, 1])[0, 1]
            if _np.abs(r2) <= max_corr:
                return values, (float(r1), float(r2))


def __modify_distributions(distr: _AbstractUnivarDistr,
                           number_list: _ntp.NDArray,
                           max_corr=0.01,) -> _Tuple[_ntp.NDArray[_np.float64], float]:

    n = len(number_list)
    while True:
        values = distr.sample(n)
        r = _np.corrcoef(number_list, values)[0, 1]
        if _np.abs(r) <= max_corr:
            return _np.atleast_2d(values).T, float(r)
