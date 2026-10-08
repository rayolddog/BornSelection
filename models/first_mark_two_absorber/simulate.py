#!/usr/bin/env python3
"""
simulate.py — continuous unitary propagation for piecewise-constant H(t).

Within one segment H is constant, so the exact solution of

    d|Psi>/dt = -(i/hbar) H |Psi>

is |Psi(t+dt)> = exp(-i H dt / hbar) |Psi(t)>.  We get that exponential from a
Hermitian eigendecomposition, which is exact for a Hermitian H up to floating
point and keeps the propagator unitary by construction.

There is no integrator tolerance to tune, no stochastic term, no jump operator
and no outcome sampling anywhere in this module.
"""
import numpy as np
from scipy.linalg import eigh, expm


def segment_propagator(H, dt, hbar=1.0):
    """exp(-i H dt / hbar) via Hermitian eigendecomposition."""
    evals, V = eigh(H)
    return (V * np.exp(-1j * evals * dt / hbar)) @ V.conj().T


def segment_propagator_expm(H, dt, hbar=1.0):
    """Same propagator via a general matrix exponential — independent cross-check."""
    return expm(-1j * H * dt / hbar)


def propagate(psi0, schedule, model, samples_per_segment=1):
    """Evolve a state through a schedule and return the sampled trajectory.

    Returns `(times, states)` where `states[0]` is the initial state and
    `states[-1]` is the final one.  Each segment is subdivided into
    `samples_per_segment` equal steps purely for output resolution; the physics
    is unaffected because H is constant inside a segment.
    """
    hbar = model.params.hbar
    psi = np.asarray(psi0, dtype=complex)
    times = [0.0]
    states = [psi.copy()]
    t = 0.0

    for seg in schedule:
        H = model.hamiltonian(seg.controls, seg.sign)
        n_steps = max(1, int(samples_per_segment))
        dt = seg.duration / n_steps
        U = segment_propagator(H, dt, hbar)
        for _ in range(n_steps):
            psi = U @ psi
            t += dt
            times.append(t)
            states.append(psi.copy())

    return np.asarray(times), np.asarray(states)


def final_state(psi0, schedule, model):
    """Just the state at the end of the schedule."""
    _, states = propagate(psi0, schedule, model)
    return states[-1]


def phase_aligned_error(psi_ref, psi):
    """Distance between two states after removing an irrelevant global phase.

    Returns `(error, overlap_magnitude)`.  The error is the 2-norm of
    `psi - exp(i phi) psi_ref` with `phi` chosen to align them, which is the
    right way to ask "did we get the complete state back?" — a global phase is
    not physical, and every other relative phase still counts.
    """
    overlap = np.vdot(psi_ref, psi)
    mag = abs(overlap)
    if mag < 1e-14:
        return float(np.linalg.norm(psi - psi_ref)), float(mag)
    aligned = (overlap / mag) * psi_ref
    return float(np.linalg.norm(psi - aligned)), float(mag)
