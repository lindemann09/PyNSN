import typing as tp
from abc import ABCMeta

from .._stimulus import NSNStimulus, NSNStimulusPair
from .._stimulus.properties import VP, VPList

ListNSNStimuli = tp.List[NSNStimulus]
ListNSNStimPairs = tp.List[NSNStimulusPair]

OptionalVPList = None | VP | str | VPList


class AbstractCollection(metaclass=ABCMeta):
    pass
