"""Exception types raised by finaleme_too."""


class FinaleMeTOOError(Exception):
    """Base exception for all finaleme_too errors."""


class InvalidSampleSheetError(FinaleMeTOOError):
    """Sample sheet is missing required columns or has invalid values."""


class InvalidInputFormatError(FinaleMeTOOError):
    """Methylation input file format cannot be parsed or detected."""


class InvalidReferencePanelError(FinaleMeTOOError):
    """Reference panel file is malformed or incompatible."""


class InvalidMarkerRegionsError(FinaleMeTOOError):
    """Marker regions file is malformed or incompatible."""


class InvalidCalibrationError(FinaleMeTOOError):
    """Calibration parameters file is malformed or missing required fields."""


class IllegalImputationError(FinaleMeTOOError):
    """Attempted to impute across comparison groups (forbidden)."""


class DeconvolutionFailedError(FinaleMeTOOError):
    """Optimization failed to converge or produced invalid output."""


class CoverageTooLowError(FinaleMeTOOError):
    """Sample coverage is below the minimum threshold for any tier."""
