LOG_SIZE = "Log Size"
TOTAL_SURFACE_AREA = "Total surface area"
ITEM_DIAMETER = "Mean item diameter"
ITEM_SURFACE_AREA = "Mean item surface area"
ITEM_PERIMETER = "Total perimeter"
TOTAL_PERIMETER = "Mean item perimeter"
LOG_SPACING = "Log Spacing"
SPARSITY = "Sparsity"
FIELD_AREA = "Field area"
COVERAGE = "Coverage"

SIZE_FEATURES = (LOG_SIZE, TOTAL_SURFACE_AREA, ITEM_DIAMETER,
                 ITEM_SURFACE_AREA, ITEM_PERIMETER, TOTAL_PERIMETER)

SPACE_FEATURES = (LOG_SPACING, SPARSITY, FIELD_AREA)

ALL_FEATURES = SIZE_FEATURES + SPACE_FEATURES + (COVERAGE,)

def are_dependent(featureA, featureB):
    """returns true if both features are not independent"""
    for l in [SIZE_FEATURES, SPACE_FEATURES]:
        if featureA in l and featureB in l:
            return True
    return False
