import typing as tp
from abc import ABCMeta

from .._stimulus import NSNStimulus, NSNStimulusPair

ListNSNStimuli = tp.List[NSNStimulus]
ListNSNStimPairs = tp.List[NSNStimulusPair]


class AbstractCollection(metaclass=ABCMeta):
    pass
