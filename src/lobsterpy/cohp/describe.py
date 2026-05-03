"""Deprecated class for describing COHPs/COOPs/COBIs."""

from monty.dev import deprecated

from lobsterpy.coxx.describe import Description as _Description


@deprecated(
    message="use `lobsterpy.coxx.describe.Description` instead.",
    deadline=(2026, 6, 30),
)
class Description(_Description):
    """Deprecated class for describing COHPs/COOPs/COBIs."""
