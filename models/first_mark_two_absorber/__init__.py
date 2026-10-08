"""
first_mark_two_absorber — the 576-state two-absorber first-mark model.

Rung 1, no-bath slice: explicit operator construction, closed unitary evolution
under a piecewise-constant Hamiltonian, the analytic capture/return control, and
the exact global inverse.

Nothing in this package samples an outcome, applies a jump, or conditions on a
result, and nothing in it supports a Born-selection or one-world claim.  See
README.md for what the model does and does not establish.
"""
from .model import (DIM, DIMS, SUBSYSTEMS, Model, Params, basis_state,
                    initial_state, kron_op, swap_AB)
from .controls import (Controls, ResonantReturnError, Segment,
                       capture_and_return_schedule, capture_schedule,
                       full_capture_time, global_inverse, mixed_schedule,
                       route_schedule, total_time)
from .simulate import (final_state, phase_aligned_error, propagate,
                       segment_propagator)

__all__ = [
    "DIM", "DIMS", "SUBSYSTEMS", "Model", "Params", "basis_state",
    "initial_state", "kron_op", "swap_AB",
    "Controls", "ResonantReturnError", "Segment", "capture_and_return_schedule",
    "capture_schedule", "full_capture_time", "global_inverse", "mixed_schedule",
    "route_schedule", "total_time",
    "final_state", "phase_aligned_error", "propagate", "segment_propagator",
]
