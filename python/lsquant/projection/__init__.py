"""Spectral window projection and orthogonal downfolding tools."""
from .spectral import (SpectralProjection, spectral_projection, resolve_projection,
                       word_trace, word_error)

__all__ = ["SpectralProjection", "spectral_projection", "resolve_projection",
           "word_trace", "word_error"]
from .orthogonal import downfold_orthogonal, validate_bands, run_origin_downfolder
from .real_space import LocalMatrix

__all__ += ["downfold_orthogonal", "validate_bands", "run_origin_downfolder", "LocalMatrix"]
