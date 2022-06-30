
DEFAULT_SPACING_PRECISION = 0.0001
DEFAULT_ADAPT_FA2TA_RATIO = 0.5

def change_adapt_settings(default_spacing_precision=None,
                          default_adapt_fa2ta_ratio=None):
    """Changing class settings of feature adapter.

    This changes the settings of all feature adapter.


    Parameters
    ----------
    default_spacing_precision
    default_adapt_fa2ta_ratio

    Returns
    -------

    """
    global DEFAULT_ADAPT_FA2TA_RATIO
    global DEFAULT_SPACING_PRECISION
    if isinstance(default_spacing_precision, float):
        DEFAULT_SPACING_PRECISION = default_spacing_precision
    if isinstance(default_adapt_fa2ta_ratio, float):
        DEFAULT_ADAPT_FA2TA_RATIO = default_adapt_fa2ta_ratio


