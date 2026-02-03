"""Deprecated class for describing COHPs/COOPs/COBIs."""

from monty.dev import deprecated

from lobsterpy.coxx.describe import Description as _Description


@deprecated(
    message="`lobsterpy.cohp.describe.Description` is deprecated; use `lobsterpy.coxx.describe.Description` instead.",
    deadline=(2026, 3, 31),
)
class Description(_Description):
    """Deprecated class for describing COHPs/COOPs/COBIs."""
