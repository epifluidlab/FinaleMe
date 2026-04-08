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
    """Calibration parameters file is malformed or missing required fields.

    DEPRECATED in v3: retained as a base class for ``InvalidBinarizationError``
    and for any v2 calibration code paths still alive during migration. New
    code should raise ``InvalidBinarizationError`` or
    ``InvalidMatchedDataError``.
    """


class InvalidBinarizationError(InvalidCalibrationError):
    """Binarization parameters file is malformed or missing required fields.

    Subclasses ``InvalidCalibrationError`` so existing ``except InvalidCalibrationError``
    blocks still catch v3 errors during the migration.
    """


class InvalidMatchedDataError(InvalidCalibrationError):
    """Matched WGBS / FinaleMe training table is malformed or unparseable.

    Shared between v2 calibration training and v3 binarization training.
    Subclasses ``InvalidCalibrationError`` so existing ``except`` blocks in v2
    code keep working during the migration.
    """


class IllegalImputationError(FinaleMeTOOError):
    """Attempted to impute across comparison groups (forbidden)."""


class DeconvolutionFailedError(FinaleMeTOOError):
    """Optimization failed to converge or produced invalid output."""


class CoverageTooLowError(FinaleMeTOOError):
    """Sample coverage is below the minimum threshold for any tier."""
