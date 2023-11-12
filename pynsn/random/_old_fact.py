# calculates visual properties of a nsn stimulus/ dot cloud

from collections import OrderedDict
from typing import Any, Optional, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _stimulus, constants
from .._lib import geometry, rng
from .._lib.exceptions import NoSolutionError
from ._object_arrays.convex_hull import ConvexHull, ConvexHullPositions
from .pynsn.constants import VisualPropertyFlags



class NSNFactory(object):



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

    def add_random_dots(self, nsn_stimulus: NSNStimulus,
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
        if not isinstance(nsn_stimulus.objects, DotArray):
            raise ValueError("nsn stimulus has to have a DotArray and not "
                             f" {type(nsn_stimulus.objects).__name__}")
        for dot in self._sample_dots(n=n_objects):
            try:
                dot = nsn_stimulus.get_free_position(ref_object=dot, in_neighborhood=False,
                                                     occupied_space=occupied_space,
                                                     allow_overlapping=allow_overlapping)
            except NoSolutionError as err:
                raise NoSolutionError(f"Can't find a solution for {n_objects} "
                                      + "items in this array") from err
            nsn_stimulus.objects.add([dot])  # type: ignore

        return nsn_stimulus

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
        rtn = NSNStimulus(object_array=DotArray(),
                          target_area_shape=deepcopy(self.target_area_shape),
                          min_dist_between_objects=self.min_dist_between_objects,
                          min_dist_area_edge=self.min_dist_area_edge)
        return self.add_random_dots(nsn_stimulus=rtn,
                                    n_objects=n_objects,
                                    allow_overlapping=allow_overlapping,
                                    occupied_space=occupied_space)

    def add_random_rectangles(self,
                              nsn_stimulus: NSNStimulus,
                              n_objects: int = 1,
                              allow_overlapping: bool = False,
                              occupied_space: Union[None, NSNStimulus] = None) -> NSNStimulus:

        if not isinstance(nsn_stimulus.objects, RectangleArray):
            raise ValueError("nsn stimulus has to have a RectangleArray and not "
                             f" {type(nsn_stimulus.objects).__name__}")
        for rect in self._sample_rectangles(n=n_objects):
            try:
                rect = nsn_stimulus.get_free_position(ref_object=rect, in_neighborhood=False,
                                                      occupied_space=occupied_space,
                                                      allow_overlapping=allow_overlapping)
            except NoSolutionError as err:
                raise NoSolutionError(f"Can't find a solution for {n_objects} "
                                      + "items in this array.") from err
            nsn_stimulus.objects.add([rect])  # type:ignore

        return nsn_stimulus

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
        rtn = NSNStimulus(object_array=RectangleArray(),
                          target_area_shape=deepcopy(self.target_area_shape),
                          min_dist_between_objects=self.min_dist_between_objects,
                          min_dist_area_edge=self.min_dist_area_edge)
        return self.add_random_rectangles(nsn_stimulus=rtn,
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
        return dict_to_text(self.to_dict())



class ArrayProperties(object):
    """Visual properties fo an associated arrays

    Visual features of each nsn stimulus can be access and
    modified via an instance of this class
    """

    def __init__(
        self, object_array: Union[ObjectArrayType, _stimulus.NSNStimulus]
    ) -> None:
        if isinstance(object_array, _stimulus.NSNStimulus):
            self._oa = object_array.objects
            self._nsn_stimulus = object_array
        elif isinstance(object_array, (DotArray, RectangleArray)):
            self._oa = object_array
            self._nsn_stimulus = None
        else:
            raise TypeError(
                "object_array has to be a NSNStimulus or an ObjectArrayType"
            )

        self._convex_hull = None
        self._convex_hull_positions = None

    def fit_numerosity(self, value: int, keep_convex_hull: bool = False) -> None:
        """ """
        if self._nsn_stimulus is None:
            raise RuntimeError(
                "Numerosity can be only modified, if NSNStimulus is "
                "defined the property, because free position can not be  "
                "found otherwise."
            )

        # make a copy for the deviant
        if value <= 0:
            self._oa.clear()
        else:
            # add or remove random dots
            change_numerosity = value - self.numerosity
            for _ in range(abs(change_numerosity)):
                if change_numerosity < 0:
                    # remove dots
                    if keep_convex_hull:
                        # find a random object that is not in convex hull
                        delete_id = None
                        ch = self.convex_hull.object_indices
                        rnd_seq = list(range(0, self.numerosity))
                        rng.generator.shuffle(rnd_seq)
                        for x in rnd_seq:
                            if x not in ch:
                                delete_id = x
                                break
                        if delete_id is None:
                            raise NoSolutionError(
                                "Can't increase numerosity, while keeping field area."
                            )
                    else:
                        delete_id = rng.generator.integers(0, self.numerosity)

                    self._oa.delete(delete_id)

                else:
                    # add object: copy a random dot
                    clone_id = rng.generator.integers(0, self.numerosity)
                    rnd_object = next(self._oa.iter(clone_id))
                    try:
                        rnd_object = self._nsn_stimulus.get_free_position(
                            ref_object=rnd_object,
                            allow_overlapping=False,
                            inside_convex_hull=keep_convex_hull,
                        )
                    except NoSolutionError as err:
                        # no free position
                        raise NoSolutionError(
                            "Can't increase numerosity. No free position found."
                        ) from err

                    self._oa.add([rnd_object])  # type: ignore

    def fit_average_diameter(self, value: float) -> None:
        """Set average diameter.

        Args:
            value: diameter

        Raises:
            TypeError: if associated array is not a array of dots
        """
        # changes diameter
        if not isinstance(self._oa, DotArray):
            raise TypeError(
                "Adapting diameter is not possible " +
                f"for {type(self._oa).__name__}."
            )
        scale = value / self.average_dot_diameter
        self._oa.diameter = self._oa.diameter * scale  # pylint: disable=W0212
        self.reset()

    def fit_average_rectangle_size(self, value: ArrayLike) -> None:
        """Set average rectangle size.

        Args:
            value:  (width, height)

        Raises:
            TypeError: if associated array is not a object of rectangles
            or values is not a tuple of two numerical values
        """
        if not isinstance(self._oa, RectangleArray):
            raise TypeError(
                "Adapting rectangle size is not possible "
                + f"for {type(self._oa).__name__}."
            )
        new_size = np.asarray(value)
        if new_size.shape != (2,):
            raise ValueError(
                f"Value ({value}) has to tuple of 2 numerical " +
                "(width, height)."
            )

        av_size = self.average_rectangle_size
        if not np.all(av_size > 0):
            raise RuntimeError(
                "Numerosity, width or hight is zero or not defined.")
        scale = np.divide(new_size, av_size)
        self._oa.sizes = self._oa.sizes * scale  # pylint: disable=W0212
        self.reset()

    def fit_total_surface_area(self, value: float) -> None:
        """Set surface area.

        Resize all object to fit a specific surface area

        Args:
            value: surface area
        """
        a_scale = value / self.total_surface_area
        # pylint: disable=W0212
        if isinstance(self._oa, DotArray):
            self._oa.diameter = (
                np.sqrt(self._oa.surface_areas * a_scale) * 2 / np.sqrt(np.pi)
            )  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
        elif isinstance(self._oa, RectangleArray):
            # rect
            self._oa.sizes = self._oa.sizes * np.sqrt(a_scale)
        else:
            raise NotImplementedError()

        self.reset()

    def fit_field_area(self, value: float, precision: Optional[float] = None) -> None:
        """changes the convex hull area to a desired size with certain precision

        uses scaling radial positions if field area has to be increased
        uses replacement of outer points (and later re-scaling)

        iterative method can takes some time.
        """

        if precision is None:
            precision = constants.DEFAULT_FIT_SPACING_PRECISION

        if self.field_area is None:
            return None  # not defined
        else:
            _match_field_area(self, value=value, precision=precision)

    def fit_coverage(
        self,
        value: float,
        precision: Optional[float] = None,
        fa2ta_ratio: Optional[float] = None,
    ) -> None:
        """

        Parameters
        ----------
        value
        precision
        fa2ta_ratio

        Returns
        -------

        """

        # TODO check drifting outwards if extra space is small and adapt_FA2TA_ratio=1
        # when to realign, realignment changes field_area!
        # """this function changes the area and remixes to get a desired density
        # precision in percent between 1 < 0
        #
        # ratio_area_convex_hull_adaptation:
        #    ratio of adaptation via area or via convex_hull (between 0 and 1)

        print("WARNING: _adapt_coverage is a experimental ")
        # dens = convex_hull_area / total_surface_area
        if fa2ta_ratio is None:
            fa2ta_ratio = constants.DEFAULT_FIT_FA2TA_RATIO
        elif fa2ta_ratio < 0 or fa2ta_ratio > 1:
            fa2ta_ratio = 0.5
        if precision is None:
            precision = constants.DEFAULT_FIT_SPACING_PRECISION

        total_area_change100 = (value * self.field_area) - \
            self.total_surface_area
        d_change_total_area = total_area_change100 * (1 - fa2ta_ratio)
        if abs(d_change_total_area) > 0:
            self.fit_total_surface_area(
                self.total_surface_area + d_change_total_area)

        self.fit_field_area(self.total_surface_area /
                            value, precision=precision)

    def fit_average_perimeter(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        self.fit_total_perimeter(value * self.numerosity)

    def fit_total_perimeter(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        if isinstance(self._oa, DotArray):
            self.fit_average_diameter(value / (self.numerosity * np.pi))

        elif isinstance(self._oa, RectangleArray):
            new_size = self.average_rectangle_size * value / self.total_perimeter

            self.fit_average_rectangle_size(new_size)

    def fit_average_surface_area(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        self.fit_total_surface_area(self.numerosity * value)

    def fit_log_spacing(self, value: float, precision: Optional[float] = None) -> None:
        """

        Parameters
        ----------
        value
        precision

        Returns
        -------

        """
        logfa = 0.5 * value + 0.5 * np.log2(self.numerosity)
        self.fit_field_area(value=2**logfa, precision=precision)

    def fit_log_size(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        logtsa = 0.5 * value + 0.5 * np.log2(self.numerosity)
        self.fit_total_surface_area(2**logtsa)

    def fit_sparcity(self, value: float, precision=None) -> None:
        """

        Parameters
        ----------
        value
        precision

        Returns
        -------

        """
        return self.fit_field_area(value=value * self.numerosity, precision=precision)

    def fit(self, property_flag: VisualPropertyFlags, value: float) -> Any:
        """
        adapt_properties: continuous property or list of continuous properties
        several properties to be adapted
        if adapt nsn stimulus is specified, array will be adapt to adapt_dot_array, otherwise
        the values defined in adapt_properties is used.
        some adapting requires realignement to avoid overlaps. However,
        realigment might result in a different field area. Thus, realign after
        adapting for  Size parameter and realign before adapting space
        parameter.

        """

        # type check
        if not isinstance(property_flag, VisualPropertyFlags):
            raise ValueError(
                "{} is not a visual feature.".format(property_flag))

        # Adapt
        if property_flag == VisualPropertyFlags.AV_DOT_DIAMETER:
            return self.fit_average_diameter(value=value)

        elif property_flag == VisualPropertyFlags.NUMEROSITY:
            return self.fit_numerosity(value=int(value))

        elif property_flag == VisualPropertyFlags.AV_PERIMETER:
            return self.fit_average_perimeter(value=value)

        elif property_flag == VisualPropertyFlags.TOTAL_PERIMETER:
            return self.fit_total_perimeter(value=value)

        elif property_flag == VisualPropertyFlags.AV_SURFACE_AREA:
            return self.fit_average_surface_area(value=value)

        elif property_flag == VisualPropertyFlags.TOTAL_SURFACE_AREA:
            return self.fit_total_surface_area(value=value)

        elif property_flag == VisualPropertyFlags.LOG_SIZE:
            return self.fit_log_size(value=value)

        elif property_flag == VisualPropertyFlags.LOG_SPACING:
            return self.fit_log_spacing(value=value)

        elif property_flag == VisualPropertyFlags.SPARSITY:
            return self.fit_sparcity(value=value)

        elif property_flag == VisualPropertyFlags.FIELD_AREA:
            return self.fit_field_area(value=value)

        elif property_flag == VisualPropertyFlags.COVERAGE:
            return self.fit_coverage(value=value)
        else:
            raise NotImplementedError(
                "Not implemented for {}".format(property_flag.label())
            )

    def scale_average_diameter(self, factor: float) -> None:
        if not isinstance(self._oa, DotArray):
            raise TypeError(
                "Scaling diameter is not possible for {}.".format(
                    type(self._oa).__name__
                )
            )
        if factor == 1:
            return
        return self.fit_average_diameter(self.average_dot_diameter * factor)

    def scale_average_rectangle_size(self, factor: float) -> None:
        if not isinstance(self._oa, RectangleArray):
            raise TypeError(
                "Scaling rectangle size is not possible for {}.".format(
                    type(self._oa).__name__
                )
            )
        if factor == 1:
            return
        return self.fit_average_rectangle_size(self.average_rectangle_size * factor)

    def scale_total_surface_area(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_total_surface_area(self.total_surface_area * factor)

    def scale_field_area(
        self, factor: float, precision: Optional[float] = None
    ) -> None:
        if factor == 1:
            return
        return self.fit_field_area(self.field_area * factor, precision=precision)

    def scale_coverage(
        self,
        factor: float,
        precision: Optional[float] = None,
        FA2TA_ratio: Optional[float] = None,
    ) -> None:
        if factor == 1:
            return
        return self.fit_coverage(
            self.coverage * factor, precision=precision, fa2ta_ratio=FA2TA_ratio
        )

    def scale_average_perimeter(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_average_perimeter(self.average_perimeter * factor)

    def scale_total_perimeter(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_total_perimeter(self.total_perimeter * factor)

    def scale_average_surface_area(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_average_surface_area(self.average_surface_area * factor)

    def scale_log_spacing(
        self, factor: float, precision: Optional[float] = None
    ) -> None:
        if factor == 1:
            return
        return self.fit_log_spacing(self.log_spacing * factor, precision=precision)

    def scale_log_size(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_log_size(self.log_size * factor)

    def scale_sparcity(self, factor: float, precision: Optional[float] = None) -> None:
        if factor == 1:
            return
        return self.fit_sparcity(self.sparsity * factor, precision=precision)

    def scale(self, feature: VisualPropertyFlags, factor: float) -> None:
        if factor == 1:
            return
        return self.fit(property_flag=feature, value=self.get(feature) * factor)


def _match_field_area(
    properties: ArrayProperties, value: float, precision: float
) -> None:
    """change the convex hull area to a desired size by scale the polar
    positions with certain precision

    iterative method can takes some time.

    Note: see doc string `field_area`
    """

    object_array = properties._oa
    current = properties.field_area

    if current is None:
        return  # not defined

    scale = 1  # find good scale
    step = 0.1
    if value < current:  # current too larger
        step *= -1

    # centered points
    old_center = properties.convex_hull.center
    object_array.xy = object_array.xy - old_center
    centered_polar = geometry.cartesian2polar(object_array.xy)

    # iteratively determine scale
    while abs(current - value) > precision:
        scale += step

        object_array.xy = geometry.polar2cartesian(centered_polar * [scale, 1])
        properties.reset()  # required at this point to recalc convex hull
        current = properties.field_area

        if (current < value and step < 0) or (current > value and step > 0):
            step *= -0.2  # change direction and finer grain

    object_array.xy = object_array.xy + old_center
    properties.reset()

    # TODO "visual test" (eye inspection) of fitting rect _arrays
