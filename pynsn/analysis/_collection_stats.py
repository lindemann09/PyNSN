from math import ceil
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress

from .._stimulus.properties import VP, ensure_vp, VPOrList
from ..collections._coll_stim_pairs import CollectionStimulusPairs


def _regression_plot(ax: plt.Axes, df: pd.DataFrame, x: str, y: str):
    slope, intercept, r_value, p_value, std_err = linregress(df[x], df[y])

    ax.scatter(df[x], df[y])
    # Add regression line
    reg_line = slope * df[x] + intercept  # type: ignore
    ax.plot(df[x], reg_line, color='green')
    ax.set_title(f"{x}, {y}  (r={r_value:.2f})")
    return


def property_regression_sheet(stim_set: CollectionStimulusPairs,
                              dv: Union[str, VP],
                              iv: VPOrList,
                              ratios: bool = False,
                              figsize: Tuple[float, float] = (10, 8)):
    if isinstance(iv, List):
        iv = [ensure_vp(p) for p in iv]
    else:
        iv = [ensure_vp(iv)]
    dv = ensure_vp(dv)

    if ratios:
        r = stim_set.property_ratios(iv + [dv])
    else:
        r = stim_set.property_differences(iv + [dv])

    if len(iv) == 1:
        n_col = 1
        n_row = 1
    else:
        n_col = 2
        n_row = ceil(len(iv)/2)

    fig, axs = plt.subplots(n_row, n_col, figsize=figsize)
    for i, p in enumerate(iv):
        _regression_plot(axs.flat[i], r, dv.name, p.name)

    plt.tight_layout()
    return fig
