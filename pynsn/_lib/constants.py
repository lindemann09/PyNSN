import enum as _enum
from .._lib import lib_typing as _tp

# _arrays
DEFAULT_MIN_DISTANCE_BETWEEN_OBJECTS = 2
DEFAULT_MIN_DISTANCE_AREA_BOARDER = 1

# iterations
MAX_ITERATIONS = 1000

# fitting
DEFAULT_FIT_SPACING_PRECISION = 0.0001
DEFAULT_FIT_FA2TA_RATIO = 0.5


class VisualPropertyFlags(_enum.IntFlag):
    AV_DOT_DIAMETER = _enum.auto()
    AV_SURFACE_AREA = _enum.auto()
    AV_PERIMETER = _enum.auto()
    AV_RECT_SIZE = _enum.auto()

    TOTAL_SURFACE_AREA = _enum.auto()
    TOTAL_PERIMETER = _enum.auto()
    SPARSITY = _enum.auto()
    FIELD_AREA = _enum.auto()
    FIELD_AREA_POSITIONS = _enum.auto()
    COVERAGE = _enum.auto()

    LOG_SPACING = _enum.auto()
    LOG_SIZE = _enum.auto()

    NUMEROSITY = _enum.auto()

    def is_dependent_from(self, other_property: _tp.Any) -> bool:
        """returns true if both properties are not independent"""
        return (self.is_size_property() and other_property.is_size_property()) or \
               (self.is_space_property() and other_property.is_space_property())

    def is_size_property(self) -> bool:
        return self in (VisualPropertyFlags.LOG_SIZE,
                        VisualPropertyFlags.TOTAL_SURFACE_AREA,
                        VisualPropertyFlags.AV_DOT_DIAMETER,
                        VisualPropertyFlags.AV_SURFACE_AREA,
                        VisualPropertyFlags.AV_PERIMETER,
                        VisualPropertyFlags.TOTAL_PERIMETER)

    def is_space_property(self) -> bool:
        return self in (VisualPropertyFlags.LOG_SPACING,
                        VisualPropertyFlags.SPARSITY,
                        VisualPropertyFlags.FIELD_AREA)

    def label(self) -> str:
        labels = {
            VisualPropertyFlags.NUMEROSITY: "Numerosity",
            VisualPropertyFlags.LOG_SIZE: "Log Size",
            VisualPropertyFlags.TOTAL_SURFACE_AREA: "Total surface area",
            VisualPropertyFlags.AV_DOT_DIAMETER: "Average dot diameter",
            VisualPropertyFlags.AV_SURFACE_AREA: "Average surface area",
            VisualPropertyFlags.AV_PERIMETER: "Average perimeter",
            VisualPropertyFlags.TOTAL_PERIMETER: "Total perimeter",
            VisualPropertyFlags.AV_RECT_SIZE: "Average Rectangle Size",
            VisualPropertyFlags.LOG_SPACING: "Log Spacing",
            VisualPropertyFlags.SPARSITY: "Sparsity",
            VisualPropertyFlags.FIELD_AREA: "Field area",
            VisualPropertyFlags.FIELD_AREA_POSITIONS: "Field area positions",
            VisualPropertyFlags.COVERAGE: "Coverage"}
        return labels[self]
