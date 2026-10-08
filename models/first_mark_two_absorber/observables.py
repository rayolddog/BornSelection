#!/usr/bin/env python3
"""
observables.py — populations, energies and overlaps read off the full state.

Everything here is a plain expectation value in the unconditioned state.  No
observable in this module selects, samples, or conditions on an outcome, and
none of them is a "which absorber won" variable.
"""
import numpy as np

from .model import (AXIS, DIMS, ABS_E, ABS_G, ABS_M, FIELD_A, FIELD_B,
                    FIELD_OUT, FIELD_VAC, QB_1)

_FIELD_NAMES = {FIELD_VAC: "vacuum", FIELD_A: "photon_A",
                FIELD_B: "photon_B", FIELD_OUT: "out"}
_ABS_NAMES = {ABS_G: "ground", ABS_E: "excited", ABS_M: "mark"}


def _tensor(psi):
    return np.asarray(psi, dtype=complex).reshape(DIMS)


def _marginal(psi, axis):
    """Population of each basis state of one subsystem."""
    probs = np.abs(_tensor(psi)) ** 2
    other = tuple(k for k in range(len(DIMS)) if k != axis)
    return probs.sum(axis=other)


def norm(psi):
    return float(np.linalg.norm(psi))


def field_populations(psi):
    """{'vacuum': .., 'photon_A': .., 'photon_B': .., 'out': ..}"""
    marg = _marginal(psi, AXIS["field"])
    return {name: float(marg[idx]) for idx, name in _FIELD_NAMES.items()}


def absorber_populations(psi, arm):
    """{'ground': .., 'excited': .., 'mark': ..} for absorber A or B."""
    marg = _marginal(psi, AXIS[f"absorber_{arm}"])
    return {name: float(marg[idx]) for idx, name in _ABS_NAMES.items()}


def record_excitation(psi, arm):
    """Probability that the local energy-record qubit is excited."""
    return float(_marginal(psi, AXIS[f"record_{arm}"])[QB_1])


def phase_qubit_z(psi, arm):
    """<Z> of a phase-history qubit."""
    marg = _marginal(psi, AXIS[f"phase_{arm}"])
    return float(marg[0] - marg[1])


def expectation(psi, operator):
    """Real expectation value <psi|operator|psi>."""
    psi = np.asarray(psi, dtype=complex)
    return float(np.real(np.vdot(psi, operator @ psi)))


def energy_ledger(psi, model, controls):
    """Expectation of every named term of H(t), plus the total.

    The point of listing the terms separately is that work done by a
    time-dependent control is visible as a change in the controlled terms rather
    than hiding inside a single number.
    """
    terms = model.named_terms(controls)
    ledger = {name: expectation(psi, op) for name, op in terms.items()}
    ledger["total"] = expectation(psi, model.hamiltonian(controls))
    return ledger


def excitation_number(psi, model):
    """Photon plus absorber excitation count, N.

    The 576-dimensional space does contain N = 2 basis states — `|vac, e_A, e_B>`
    and `|photon-A, e_A, g_B>` are both perfectly good vectors in it.  What keeps
    the calibration runs at N = 1 is that the prepared state starts there and
    H conserves N, not the size or truncation of the space.

    Either way this is bookkeeping.  It says nothing about whether one
    alternative rather than another is realised.
    """
    return expectation(psi, model.N_excitation)


def interference_visibility(psi):
    """Fringe visibility of the two returned arm amplitudes.

    Computed from the reduced field coherence between arm A and arm B, so it
    reports how much A/B relative phase information is still accessible in the
    field alone.

    Returns NaN — not zero — when there is no photon amplitude left in either
    arm.  While the photon is held inside the absorbers the fringe visibility of
    the arm field is simply undefined, and reporting it as zero would read as
    lost coherence when nothing has been lost.
    """
    t = _tensor(psi)
    flat = t.reshape(DIMS[0], -1)
    rho_AA = float(np.vdot(flat[FIELD_A], flat[FIELD_A]).real)
    rho_BB = float(np.vdot(flat[FIELD_B], flat[FIELD_B]).real)
    rho_AB = complex(np.vdot(flat[FIELD_A], flat[FIELD_B]))
    denom = rho_AA + rho_BB
    if denom < 1e-12:
        return float("nan")
    return float(2.0 * abs(rho_AB) / denom)


def arm_relative_phase(psi):
    """Relative phase between the arm-A and arm-B field amplitudes.

    Only meaningful while the field part is still an unentangled product with
    the rest; in the no-bath calibration that is the case, which is exactly why
    the no-bath slice can check phase-faithful return directly.
    """
    t = _tensor(psi)
    flat = t.reshape(DIMS[0], -1)
    overlap = complex(np.vdot(flat[FIELD_A], flat[FIELD_B]))
    if abs(overlap) < 1e-14:
        return float("nan")
    return float(np.angle(overlap))


def summary(psi, model, controls=None):
    """A compact dict of everything worth printing at one interruption time."""
    out = {
        "norm": norm(psi),
        "field": field_populations(psi),
        "absorber_A": absorber_populations(psi, "A"),
        "absorber_B": absorber_populations(psi, "B"),
        "record_A": record_excitation(psi, "A"),
        "record_B": record_excitation(psi, "B"),
        "phase_Z_A": phase_qubit_z(psi, "A"),
        "phase_Z_B": phase_qubit_z(psi, "B"),
        "excitation_number": excitation_number(psi, model),
        "visibility": interference_visibility(psi),
    }
    if controls is not None:
        out["energy"] = energy_ledger(psi, model, controls)
    return out
