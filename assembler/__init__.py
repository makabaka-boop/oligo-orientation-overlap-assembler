"""Short-read DNA assembler: optimal order + strand selection."""

from .core import (
    Assembly,
    AssemblyError,
    Placement,
    Read,
    assemble,
    assembly_to_data,
    max_overlap,
    reads_from_data,
    reverse_complement,
    validate_reads,
)

__all__ = [
    "Assembly",
    "AssemblyError",
    "Placement",
    "Read",
    "assemble",
    "assembly_to_data",
    "max_overlap",
    "reads_from_data",
    "reverse_complement",
    "validate_reads",
]

__version__ = "1.0.0"
