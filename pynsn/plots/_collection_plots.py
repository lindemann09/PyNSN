from math import ceil
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.stats import linregress

from .._stimulus.properties import VP, VPOrList, ensure_vp
from ..collections._coll_stim_pairs import (CollectionStimuli,
                                            CollectionStimulusPairs)
from .. import defaults


def scatter_matrix(stimuli: CollectionStimuli,
                   props: VPOrList | None = None,
                   figsize: Tuple[float, float] = (10, 8),
                   alpha: float = 0.7) -> Tuple[Figure, Axes]:

    if props is None:
        cols = defaults.PLOT_PROPERTY_LIST
    elif isinstance(props, List):
        cols = [ensure_vp(p).name for p in props]
    else:
        cols = [ensure_vp(props).name]

    df = stimuli.property_dataframe()[cols]

    fig = plt.figure(figsize=figsize)
    ax = fig.subplots(len(cols), len(cols))
    pd.plotting.scatter_matrix(df, ax=ax, alpha=alpha)
    return fig, ax


def property_regression(stimuli: CollectionStimuli,
                        dv: str | VP = "N",
                        iv: VPOrList | None = None,
                        figsize: Tuple[float, float] = (10, 8,),
                        alpha: float = 0.7) -> Tuple[Figure, Axes]:

    if iv is None:
        iv_prop = [VP[p] for p in defaults.PLOT_PROPERTY_LIST]
    elif isinstance(iv, List):
        iv_prop = [ensure_vp(p) for p in iv]
    else:
        iv_prop = [ensure_vp(iv)]
    dv = ensure_vp(dv)

    data = stimuli.property_dataframe()
    cols = [p.name for p in iv_prop + [dv]]

    return _regression_plots(data[cols], dv=dv, ivs=iv_prop,
                             figsize=figsize, alpha=alpha)


def property_difference_regression(stim_pairs: CollectionStimulusPairs,
                                   dv: str | VP = "N",
                                   iv: VPOrList | None = None,
                                   figsize: Tuple[float, float] = (10, 8),
                                   alpha: float = 0.7) -> Tuple[Figure, Axes]:

    if iv is None:
        iv_prop = [VP[p] for p in defaults.PLOT_PROPERTY_LIST]
    elif isinstance(iv, List):
        iv_prop = [ensure_vp(p) for p in iv]
    else:
        iv_prop = [ensure_vp(iv)]
    dv = ensure_vp(dv)

    data = stim_pairs.property_differences(iv_prop + [dv])
    return _regression_plots(data, dv=dv, ivs=iv_prop,
                             figsize=figsize, alpha=alpha)


def property_ratio_regression(stim_pairs: CollectionStimulusPairs,
                              dv: str | VP = "N",
                              iv: VPOrList | None = None,
                              figsize: Tuple[float, float] = (10, 8),
                              alpha: float = 0.7) -> Tuple[Figure, Axes]:

    if iv is None:
        iv_prop = [VP[p] for p in defaults.PLOT_PROPERTY_LIST]
    elif isinstance(iv, List):
        iv_prop = [ensure_vp(p) for p in iv]
    else:
        iv_prop = [ensure_vp(iv)]
    dv = ensure_vp(dv)

    data = stim_pairs.property_ratios(iv_prop + [dv])
    return _regression_plots(data, dv=dv, ivs=iv_prop,
                             figsize=figsize, alpha=alpha)


def _regression_plots(data: pd.DataFrame,
                      dv:  VP,
                      ivs: List[VP],
                      figsize: Tuple[float, float],
                      alpha: float) -> Tuple[Figure, Axes]:

    if dv in ivs:
        raise ValueError(f"Dependent variable '{dv.name}' is also"
                         " in list of independent variables.")
    if len(ivs) == 1:
        n_col = 1
        n_row = 1
    else:
        n_col = 2
        n_row = ceil(len(ivs)/2)

    fig = plt.figure(figsize=figsize)
    axs = fig.subplots(n_row, n_col)

    for i, p in enumerate(ivs):
        __reg_plot(axs.flat[i], data, dv.name, p.name, alpha=alpha)

    fig.tight_layout()

    return fig, axs


def __reg_plot(ax: Axes, df: pd.DataFrame, x: str, y: str, alpha: float):  # type: ignore
    slope, intercept, r_value, _, _ = linregress(df[x], df[y])

    ax.scatter(df[x], df[y], alpha=alpha)
    # Add regression line
    reg_line = slope * df[x] + intercept  # type: ignore
    ax.plot(df[x], reg_line, color='darkgreen', alpha=0.7)
    ax.set_title(f"{x}, {y}  (r={r_value:.2f})")
