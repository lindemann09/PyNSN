"""export shape to data tables"""

from ._stimulus.nsn_stimulus import table_dict
from ._stimulus.nsn_stimulus import NSNStimulus as _NSNStimulus

def dataframe(stim:_NSNStimulus):
    """TODO"""
    import pandas
    return pandas.DataFrame(table_dict(stim))

def arrow_table(stim:_NSNStimulus):
    """TODO"""
    import pyarrow
    return pyarrow.Table.from_pydict(table_dict(stim)) # Arrow Table
