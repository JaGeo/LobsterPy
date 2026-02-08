"""Deprecated class for analyzing COHPs/COOPs/COBIs."""

from monty.dev import deprecated

from lobsterpy.coxx.analyze import Analysis as _Analysis


@deprecated(message="use `lobsterpy.coxx.analyze.Analysis` instead.", deadline=(2026, 3, 31))
class Analysis(_Analysis):
    """Deprecated class for analyzing COHPs/COOPs/COBIs."""
