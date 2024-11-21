from math import ceil
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress

from .._stimulus.properties import VP, ensure_vp, VPOrList
from ..collections._coll_stim_pairs import CollectionStimulusPairs, CollectionStimuli


def property_regression(stimuli: CollectionStimuli,
                        dv: Union[str, VP],
                        iv: VPOrList,
                        figsize: Tuple[float, float] = (10, 8)):

    if isinstance(iv, List):
        iv_prop = [ensure_vp(p) for p in iv]
    else:
        iv_prop = [ensure_vp(iv)]
    dv = ensure_vp(dv)

    data = stimuli.property_dataframe()
    cols = [p.name for p in iv_prop + [dv]]

    return _regression_plots(data[cols], dv=dv,
                             ivs=iv_prop, figsize=figsize)


def property_difference_regression(stim_pairs: CollectionStimulusPairs,
                                   dv: Union[str, VP],
                                   iv: VPOrList,
                                   figsize: Tuple[float, float] = (10, 8)):

    if isinstance(iv, List):
        iv_prop = [ensure_vp(p) for p in iv]
    else:
        iv_prop = [ensure_vp(iv)]
    dv = ensure_vp(dv)

    data = stim_pairs.property_differences(iv_prop + [dv])
    return _regression_plots(data, dv=dv,
                             ivs=iv_prop, figsize=figsize)


def property_ratio_regression(stim_pairs: CollectionStimulusPairs,
                              dv: Union[str, VP],
                              iv: VPOrList,
                              figsize: Tuple[float, float] = (10, 8)):

    if isinstance(iv, List):
        iv_prop = [ensure_vp(p) for p in iv]
    else:
        iv_prop = [ensure_vp(iv)]
    dv = ensure_vp(dv)

    data = stim_pairs.property_ratios(iv_prop + [dv])
    return _regression_plots(data, dv=dv,
                             ivs=iv_prop, figsize=figsize)


def _regression_plots(data: pd.DataFrame,
                      dv:  VP,
                      ivs: List[VP],
                      figsize: Tuple[float, float] = (10, 8)):

    if dv in ivs:
        raise ValueError(f"Dependent variable '{dv.name}' is also"
                         " in list of independent variables.")
    if len(ivs) == 1:
        n_col = 1
        n_row = 1
    else:
        n_col = 2
        n_row = ceil(len(ivs)/2)

    fig, axs = plt.subplots(n_row, n_col, figsize=figsize)
    for i, p in enumerate(ivs):
        _reg_plot(axs.flat[i], data, dv.name, p.name)

    plt.tight_layout()
    return fig


def _reg_plot(ax: plt.Axes, df: pd.DataFrame, x: str, y: str):  # type: ignore
    slope, intercept, r_value, _, _ = linregress(df[x], df[y])

    ax.scatter(df[x], df[y])
    # Add regression line
    reg_line = slope * df[x] + intercept  # type: ignore
    ax.plot(df[x], reg_line, color='green')
    ax.set_title(f"{x}, {y}  (r={r_value:.2f})")
    return
