import enum as _enum
from .._lib import lib_typing as _tp

# _arrays
DEFAULT_MIN_DIST_BETWEEN = 2
DEFAULT_MIN_DIST_AREA_BOARDER = 1

# iterations
MAX_ITERATIONS = 1000

# fitting
DEFAULT_FIT_SPACING_PRECISION = 0.0001
DEFAULT_FIT_FA2TA_RATIO = 0.5


class VisualPropertyFlag(_enum.IntFlag):
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
        return self in (VisualPropertyFlag.LOG_SIZE,
                        VisualPropertyFlag.TOTAL_SURFACE_AREA,
                        VisualPropertyFlag.AV_DOT_DIAMETER,
                        VisualPropertyFlag.AV_SURFACE_AREA,
                        VisualPropertyFlag.AV_PERIMETER,
                        VisualPropertyFlag.TOTAL_PERIMETER)

    def is_space_property(self) -> bool:
        return self in (VisualPropertyFlag.LOG_SPACING,
                        VisualPropertyFlag.SPARSITY,
                        VisualPropertyFlag.FIELD_AREA)

    def label(self) -> str:
        labels = {
            VisualPropertyFlag.NUMEROSITY: "Numerosity",
            VisualPropertyFlag.LOG_SIZE: "Log Size",
            VisualPropertyFlag.TOTAL_SURFACE_AREA: "Total surface area",
            VisualPropertyFlag.AV_DOT_DIAMETER: "Average dot diameter",
            VisualPropertyFlag.AV_SURFACE_AREA: "Average surface area",
            VisualPropertyFlag.AV_PERIMETER: "Average perimeter",
            VisualPropertyFlag.TOTAL_PERIMETER: "Total perimeter",
            VisualPropertyFlag.AV_RECT_SIZE: "Average Rectangle Size",
            VisualPropertyFlag.LOG_SPACING: "Log Spacing",
            VisualPropertyFlag.SPARSITY: "Sparsity",
            VisualPropertyFlag.FIELD_AREA: "Field area",
            VisualPropertyFlag.FIELD_AREA_POSITIONS: "Field area positions",
            VisualPropertyFlag.COVERAGE: "Coverage"}
        return labels[self]
