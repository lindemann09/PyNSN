
ITERATIVE_CONVEX_HULL_MODIFICATION = False
TAKE_RANDOM_DOT_FROM_CONVEXHULL = True  # TODO maybe method for modification, set via GUI
DEFAULT_SPACING_PRECISION = 0.0001
DEFAULT_ADAPT_FA2TA_RATIO = 0.5

def change_adapt_settings(iterative_convex_hull_modification=None,
                          take_random_dot_from_convexhull=None,
                          default_spacing_precision=None,
                          default_adapt_fa2ta_ratio=None):
    """Changing class settings of feature adapter.

    This changes the settings of all feature adapter.


    Parameters
    ----------
    iterative_convex_hull_modification
    take_random_dot_from_convexhull
    default_spacing_precision
    default_adapt_fa2ta_ratio

    Returns
    -------

    """
    global ITERATIVE_CONVEX_HULL_MODIFICATION
    global TAKE_RANDOM_DOT_FROM_CONVEXHULL
    global DEFAULT_ADAPT_FA2TA_RATIO
    global DEFAULT_SPACING_PRECISION
    if isinstance(iterative_convex_hull_modification, bool):
        ITERATIVE_CONVEX_HULL_MODIFICATION = iterative_convex_hull_modification
    if isinstance(take_random_dot_from_convexhull, bool):
        TAKE_RANDOM_DOT_FROM_CONVEXHULL = take_random_dot_from_convexhull
    if isinstance(default_spacing_precision, float):
        DEFAULT_SPACING_PRECISION = default_spacing_precision
    if isinstance(default_adapt_fa2ta_ratio, float):
        DEFAULT_ADAPT_FA2TA_RATIO = default_adapt_fa2ta_ratio


