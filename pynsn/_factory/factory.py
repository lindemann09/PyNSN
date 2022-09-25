__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import copy, deepcopy
from typing import Iterator, Optional, Sequence, Union

from .._arrays.nsn_stimulus import NSNStimulus
from .._arrays.target_area import TargetArea
from .._lib.exception import NoSolutionError
from .._lib.misc import dict_to_text
from .._shapes.dot_list import DotList
from .._shapes.rectangle_list import RectangleList
from .._shapes.dot import Dot
from .._shapes.rectangle import Rectangle
from .._shapes import ShapeType
from .distributions import Levels, PyNSNDistribution, _Constant, ParameterDistributionType


def _make_distr(value: ParameterDistributionType) -> Union[PyNSNDistribution, None]:
    """helper
    returns a distribution or None, if None
    """

    if value is None:
        return None
    elif isinstance(value, PyNSNDistribution):
        return value
    elif isinstance(value, (list, tuple)):
        return Levels(levels=copy(value))
    elif isinstance(value, (float, int)):
        return _Constant(value)
    else:
        raise RuntimeError("Can't make distribution from"
                           f" {type(value).__name__}")


class NSNFactory(TargetArea):

    def __init__(self,
                 target_area: ShapeType,
                 min_dist_between_objects: Optional[float] = None,
                 min_dist_area_edge: Optional[float] = None):
        """

        Parameters
        ----------
        target_area_radius
        min_dist_between_objects
        min_dist_area_edge
        """

        TargetArea.__init__(self,
                            target_area_shape=target_area,
                            min_dist_between_objects=min_dist_between_objects,
                            min_dist_area_edge=min_dist_area_edge)
        self._distr_diameter = None
        self._distr_width = None
        self._distr_height = None
        self._distr_proportion = None
        self._distr_dot_attributes = None
        self._distr_rect_attributes = None

    @property
    def distr_diameter(self) -> Optional[PyNSNDistribution]:
        """Distribution of diameter parameter"""
        return self._distr_diameter

    @property
    def distr_width(self) -> Optional[PyNSNDistribution]:
        """Distribution of width parameter"""
        return self._distr_width

    @property
    def distr_height(self) -> Optional[PyNSNDistribution]:
        """Distribution of height parameter"""
        return self._distr_height

    @property
    def distr_proportion(self) -> Optional[PyNSNDistribution]:
        """Distribution of proportion parameter"""
        return self._distr_proportion

    @property
    def distr_dot_attributes(self) -> Optional[PyNSNDistribution]:
        """Distribution of attributes for dots"""
        return self._distr_dot_attributes

    @property
    def distr_rectangle_attributes(self) -> Optional[PyNSNDistribution]:
        """Distribution of attributes for rectangles"""
        return self._distr_dot_attributes

    def set_appearance_dots(self,
                            diameter: ParameterDistributionType,
                            attributes=None):
        """Set distributions of the parameter for random nsn stimuli

        Args:
            diameter: distribution of diameter
            attributes: distribution of attributes (Optional)
        """
        self._distr_diameter = _make_distr(diameter)
        self._distr_dot_attributes = _make_distr(attributes)

    def set_appearance_rectangles(self,
                                  width: ParameterDistributionType = None,
                                  height: ParameterDistributionType = None,
                                  proportion: ParameterDistributionType = None,
                                  attributes: ParameterDistributionType = None):
        """Set distributions of the parameter for random nsn stimuli

        Args:
            width: distribution of width (Optional)
            height: distribution of height (Optional)
            proportion: distribution of proportions (Optional)
            attributes: distribution of attributes (Optional)

        Notes:
            Define either rectangle width and height or rectangle proportion together
            with either width or height.

        Raises:
            TypeError: if not two of the three rectangle parameter are defined
        """
        n_rect_parameter = sum([width is not None, height is not None,
                                proportion is not None])
        if n_rect_parameter != 2:
            raise TypeError("Define rectangle width and height or, alternatively, rectangle proportion together with "
                            "either width or height.")
        self._distr_width = _make_distr(width)
        self._distr_height = _make_distr(height)
        self._distr_proportion = _make_distr(proportion)
        self._distr_rect_attributes = _make_distr(attributes)

    def _sample_dots(self, n: int) -> Sequence[Dot]:
        """return list of dots with random size all positions = (0,0)"""
        if self._distr_diameter is None:
            raise RuntimeError(
                "Dot appearance is not defined; please use set_dot_appearance.")
        diameter = self._distr_diameter.sample(n)
        if self._distr_dot_attributes is not None:
            attributes = self._distr_dot_attributes.sample(n)
        else:
            attributes = [None] * n

        return [Dot(xy=(0, 0), diameter=dia, attribute=attr)
                for dia, attr in zip(diameter, attributes)]  # type: ignore

    def _sample_rectangles(self, n: int) -> Sequence[Rectangle]:
        """return list of rectangles with random size all positions = (0,0)"""
        if self._distr_width is None and self._distr_height is None:
            raise RuntimeError(
                "Rectangle appearance is not  defined; please use set_rectangle_appearance.")
        if self._distr_width is not None:
            width = self._distr_width.sample(n)
        else:
            width = None
        if self._distr_height is not None:
            height = self._distr_height.sample(n)
        else:
            height = None
        if self._distr_proportion is not None:
            proportion = self._distr_proportion.sample(n)
        else:
            proportion = None

        if height is None:
            height = width * proportion  # type: ignore
        elif width is None:
            width = height / proportion  # type: ignore

        if self._distr_dot_attributes is not None:
            attributes = self._distr_dot_attributes.sample(n)
        else:
            attributes = [None] * n

        return [Rectangle(xy=(0, 0), size=(w, h), attribute=attr)
                for w, h, attr in zip(width, height, attributes)]  # type: ignore

    def to_dict(self) -> dict:
        """Dict representation of the NSNFactory instance"""
        rtn = TargetArea.to_dict(self)
        try:
            rtn.update({"diameter":
                        self._distr_diameter.to_dict()})  # type: ignore
        except AttributeError:
            pass
        try:
            rtn.update({"width": self._distr_width.to_dict()})  # type: ignore
        except AttributeError:
            pass
        try:

            rtn.update({"height":
                        self._distr_height.to_dict()})  # type: ignore
        except AttributeError:
            pass
        try:
            rtn.update({"proportion":
                        self._distr_proportion.to_dict()})  # type: ignore
        except AttributeError:
            pass
        try:
            rtn.update({"attributes":
                        self._distr_attributes.to_dict()})  # type: ignore
        except AttributeError:
            pass
        return rtn

    def add_random_dots(self, object_array: NSNStimulus,
                        n_objects: int = 1,
                        allow_overlapping: bool = False,
                        occupied_space: Union[None, NSNStimulus] = None) -> NSNStimulus:
        """
        occupied_space is a nsn stimulus (used for multicolour nsn stimulus (join after)

        attribute is an array, _arrays are assigned randomly.


        Parameters
        ----------
        n_objects
        allow_overlapping
        occupied_space

        Returns
        -------
        rtn : nsn stimulus
        """
        if not isinstance(object_array.objects, DotList):
            raise ValueError("nsn stimulus has to have a DotList and not "
                             f" {type(object_array.objects).__name__}")
        for dot in self._sample_dots(n=n_objects):
            try:
                dot = object_array.get_free_position(ref_object=dot, in_neighborhood=False,
                                                     occupied_space=occupied_space,
                                                     allow_overlapping=allow_overlapping)
            except NoSolutionError as err:
                raise NoSolutionError(f"Can't find a solution for {n_objects} "
                                      + "items in this array") from err
            object_array.objects.add([dot])  # type: ignore

        return object_array

    def random_dot_array(self, n_objects: int,
                         allow_overlapping: bool = False,
                         occupied_space: Union[None, NSNStimulus] = None) -> NSNStimulus:
        """Create a new random nsn stimulus
        occupied_space is a nsn stimulus (used for multicolour nsn stimulus (join after)

        attribute is an array, _arrays are assigned randomly.


        Parameters
        ----------
        n_objects
        allow_overlapping
        occupied_space

        Returns
        -------
        rtn : nsn stimulus
        """
        rtn = NSNStimulus(object_list=DotList(),
                          target_area_shape=deepcopy(self.target_area_shape),
                          min_dist_between_objects=self.min_dist_between_objects,
                          min_dist_area_edge=self.min_dist_area_edge)
        return self.add_random_dots(object_array=rtn,
                                    n_objects=n_objects,
                                    allow_overlapping=allow_overlapping,
                                    occupied_space=occupied_space)

    def add_random_rectangles(self,
                              object_array: NSNStimulus,
                              n_objects: int = 1,
                              allow_overlapping: bool = False,
                              occupied_space: Union[None, NSNStimulus] = None) -> NSNStimulus:

        if not isinstance(object_array.objects, RectangleList):
            raise ValueError("nsn stimulus has to have a RectangleList and not "
                             f" {type(object_array.objects).__name__}")
        for rect in self._sample_rectangles(n=n_objects):
            try:
                rect = object_array.get_free_position(ref_object=rect, in_neighborhood=False,
                                                      occupied_space=occupied_space,
                                                      allow_overlapping=allow_overlapping)
            except NoSolutionError as err:
                raise NoSolutionError(f"Can't find a solution for {n_objects} "
                                      + "items in this array.") from err
            object_array.objects.add([rect])  # type:ignore

        return object_array

    def random_rectangle_array(self, n_objects: int,
                               allow_overlapping: bool = False,
                               occupied_space: Union[None, NSNStimulus] = None) -> NSNStimulus:
        """
        occupied_space is a nsn stimulus (used for multicolour nsn stimulus (join after)

        attribute is an array, _arrays are assigned randomly.


        Parameters
        ----------
        n_objects
        allow_overlapping
        occupied_space

        Returns
        -------
        rtn : nsn stimulus
        """
        rtn = NSNStimulus(object_list=RectangleList(),
                          target_area_shape=deepcopy(self.target_area_shape),
                          min_dist_between_objects=self.min_dist_between_objects,
                          min_dist_area_edge=self.min_dist_area_edge)
        return self.add_random_rectangles(object_array=rtn,
                                          n_objects=n_objects,
                                          allow_overlapping=allow_overlapping,
                                          occupied_space=occupied_space)

    def incremental_random_dot_array(self, n_objects: int,
                                     allow_overlapping: bool = False) -> Iterator[NSNStimulus]:
        """

        Parameters
        ----------
        n_objects
        allow_overlapping

        Returns
        -------
        rtn : iterator of object _arrays
        """
        previous = None
        for _ in range(n_objects):
            current = self.random_dot_array(n_objects=1,
                                            allow_overlapping=allow_overlapping,
                                            occupied_space=previous)
            if previous is not None:
                current.join(previous)
            previous = current
            yield current

    def incremental_random_rectangle_array(self, n_objects: int,
                                           allow_overlapping: bool = False) -> Iterator[NSNStimulus]:
        """

        Parameters
        ----------
        n_objects
        allow_overlapping

        Returns
        -------
        rtn : iterator of object _arrays
        """
        previous = None
        for n in range(n_objects):
            current = self.random_rectangle_array(n_objects=1,
                                                  allow_overlapping=allow_overlapping,
                                                  occupied_space=previous)
            if previous is not None:
                current.join(previous)
            previous = current
            yield current

    def __str__(self) -> str:
        return dict_to_text(self.to_dict(), col_a=12, col_b=7)
