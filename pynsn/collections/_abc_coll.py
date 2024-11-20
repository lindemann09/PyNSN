import typing as tp
from abc import ABCMeta

from .._stimulus import NSNStimulus, NSNStimulusPair
from .._stimulus.properties import VisProp, VisPropList

ListNSNStimuli = tp.List[NSNStimulus]
ListNSNStimPairs = tp.List[NSNStimulusPair]

OptionalVisPropList = None | VisProp | str | VisPropList


class AbstractCollection(metaclass=ABCMeta):
    pass
