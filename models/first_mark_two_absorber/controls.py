#!/usr/bin/env python3
"""
controls.py — piecewise-constant control schedules and the exact global inverse.

A `Schedule` is an ordered list of `Segment`s.  Within a segment every control is
constant, so H is constant there and the propagator is a single matrix
exponential.  A time-dependent control is represented by using more segments.

The exact global inverse is defined exactly as the planning artifact states it:

    H_global_reverse(s) = -H_forward(tau - s)

which for a piecewise-constant schedule means "same segments, reverse order,
every Hamiltonian negated".  `sign` on a segment carries that minus, and it is
applied to the whole assembled Hamiltonian — so the mark coupling, the phase
history, the record-qubit free energy and every field/absorber term all reverse
together.  That is deliberately *not* the operational reversal generator, which
is not implemented in this no-bath slice.
"""
from dataclasses import dataclass, replace

import numpy as np


@dataclass(frozen=True)
class Controls:
    """The externally set knobs, constant over one segment.

    K_i        capture coupling in arm i (forward, off, or reversed sign)
    delta_i    controlled absorber detuning
    u_i        prepared reference projection for the phase-history term,
               a real number in [-1, 1]; multiplied by chi_i from `Params`
    J_i        recombiner strength routing arm i into the output mode
    theta_i    recombination phase in arm i
    """
    K_A: float = 0.0
    K_B: float = 0.0
    delta_A: float = 0.0
    delta_B: float = 0.0
    u_A: float = 0.0
    u_B: float = 0.0
    J_A: float = 0.0
    J_B: float = 0.0
    theta_A: float = 0.0
    theta_B: float = 0.0

    def K(self, arm):
        return self.K_A if arm == "A" else self.K_B

    def delta(self, arm):
        return self.delta_A if arm == "A" else self.delta_B

    def u(self, arm):
        return self.u_A if arm == "A" else self.u_B

    def J(self, arm):
        return self.J_A if arm == "A" else self.J_B

    def theta(self, arm):
        return self.theta_A if arm == "A" else self.theta_B

    def with_capture(self, K):
        """Same controls with a symmetric capture coupling."""
        return replace(self, K_A=K, K_B=K)


@dataclass(frozen=True)
class Segment:
    """One interval of constant control.

    `sign` is +1 for forward evolution and -1 for the exact global inverse.
    """
    duration: float
    controls: Controls
    sign: float = 1.0
    label: str = ""


def total_time(schedule):
    return float(sum(seg.duration for seg in schedule))


def global_inverse(schedule):
    """The exact Loschmidt reversal: replay backward with the sign reversed.

    Nothing is left out.  This is the mathematical inverse, not an
    experimentally available control sequence.
    """
    return [Segment(seg.duration, seg.controls, -seg.sign,
                    label=f"inverse:{seg.label}" if seg.label else "inverse")
            for seg in reversed(schedule)]


# ---------------------------------------------------------------------------
# Analytic no-bath capture / return control
# ---------------------------------------------------------------------------
def full_capture_time(K):
    """Pulse length that transfers all photon amplitude into the absorbers.

    On resonance the pair (photon in arm i, absorber i excited) is an exactly
    degenerate two-level system with coupling hbar*K, so the transferred
    population is sin^2(K t) and a pi/2 pulse area (K t = pi/2) is complete
    capture.
    """
    return 0.5 * np.pi / K


class ResonantReturnError(ValueError):
    """Raised when the resonant return helper is asked for something it is not.

    Flipping the sign of K inverts the capture term and nothing else.  Any other
    term that is still running — a detuning, a phase-history coupling, a mark
    coupling — keeps evolving forward through all three segments and is simply
    not undone, so the photon does not come back.  Rather than quietly returning
    a schedule that no longer does what its name says, the helper refuses.
    """


def _resolve_params(context):
    """Accept a `Params` or a `Model` and return the `Params`.

    The physical context is mandatory rather than optional, because the mark
    coupling `G` lives in `Params` and not in `Controls`: no value a caller
    passes for the controls can reveal it.  A helper that promises a no-bath
    return therefore cannot validate that promise without being handed the model
    it will be evolved under.
    """
    params = getattr(context, "params", context)          # a Model, or a Params
    missing = [name for name in ("G_A", "G_B", "chi_A", "chi_B")
               if not hasattr(params, name)]
    if missing:
        raise TypeError(
            "capture_and_return_schedule() needs the physical context as its "
            "first argument: a Params or a Model, so that the no-bath contract "
            f"can actually be checked. Got {type(context).__name__!r}, which is "
            f"missing {', '.join(missing)}. Call it as "
            "capture_and_return_schedule(params, K, ...).")
    return params


def capture_and_return_schedule(context, K, t_dark=0.0, capture_time=None,
                                delta=0.0, u=0.0):
    """Capture the photon, hold, then return it by reversing the capture sign.

    `context` is the `Params` (or the `Model`) the schedule will be evolved
    under, and it is **required**.  See below for why it is not optional.

    Three segments: `+K` for a pi/2 pulse area, `K = 0` for an optional dark
    hold, then `-K` for the same pulse area.  Reversing the sign of K is one of
    the permitted operational controls; on resonance the return is coherent, so
    the photon comes back into its original arms with its original relative phase
    (up to an overall phase).

    **This helper is resonant and no-bath only, and it enforces that.**  The pulse
    length `pi/(2K)` and the sign flip on K together invert the *resonant*
    two-level capture evolution.  They do not invert:

    - a controlled detuning `Delta`, whose `P_e` term keeps running forward
      through capture, hold and return;
    - a phase-history coupling `chi * u`, for the same reason;
    - the mark coupling `G`, which the plan forbids operational reversal from
      touching at all.

    `delta` and `u` are accepted only as zero.  The first two of those live in
    `Controls`, so a caller reveals them just by calling.  `G` does not: it lives
    in `Params`, so a schedule built in ignorance of the model and then evolved
    with a non-zero `G` would fail to return with nothing raised anywhere.  That
    is why `context` is mandatory and why the `G` check is unconditional.  Every
    violation raises `ResonantReturnError` before any evolution happens.

    Concretely, these are the misuses being refused, at K = 0.85, t_dark = 1.9:

    | violation      | arm-photon population returned |
    | -------------- | ------------------------------ |
    | delta = 0.5    | 0.823                          |
    | G = 0.1        | 0.830                          |
    | G = 0.3        | 0.087                          |
    | G = 0.5        | 0.254                          |

    A non-zero `chi` is *not* rejected, and that is deliberate rather than an
    oversight: the phase term is `hbar * (chi * u / 2) * P_e * Z`, so forcing
    `u = 0` makes it vanish identically whatever `chi` is.  The suite verifies
    that a `chi != 0`, `G = 0` model still returns the photon exactly.

    Detuning belongs on `capture_schedule`, which makes no return claim and is
    used for the generalized-Rabi calibration.  A protocol that genuinely returns
    the photon in the presence of these terms is either the exact global inverse
    (`global_inverse`, which reverses everything) or the operational-reversal
    generator, which is not implemented in this slice.
    """
    params = _resolve_params(context)

    if params.G_A != 0.0 or params.G_B != 0.0:
        raise ResonantReturnError(
            f"capture_and_return_schedule is no-bath only, but the model has "
            f"G_A={params.G_A!r}, G_B={params.G_B!r}. The mark coupling keeps "
            "running through all three segments and operational reversal is not "
            "permitted to reverse it; only global_inverse() does.")
    if delta != 0.0:
        raise ResonantReturnError(
            f"capture_and_return_schedule is resonant only, but delta={delta!r}. "
            "Flipping the sign of K does not invert the detuning term, so the "
            "photon would not fully return. Use capture_schedule(delta=...) for "
            "detuned calibration, or global_inverse() for exact reversal.")
    if u != 0.0:
        raise ResonantReturnError(
            f"capture_and_return_schedule is no-bath only, but u={u!r}. "
            "Flipping the sign of K does not invert the phase-history term.")

    on = Controls(K_A=K, K_B=K)
    off = Controls()
    back = Controls(K_A=-K, K_B=-K)

    t_pulse = full_capture_time(K) if capture_time is None else capture_time
    schedule = [Segment(t_pulse, on, label="capture")]
    if t_dark > 0.0:
        schedule.append(Segment(t_dark, off, label="hold"))
    schedule.append(Segment(t_pulse, back, label="return"))
    return schedule


def capture_schedule(K, duration, delta=0.0, u=0.0):
    """A single capture pulse of the given length — the Rabi calibration run."""
    on = Controls(K_A=K, K_B=K, delta_A=delta, delta_B=delta, u_A=u, u_B=u)
    return [Segment(duration, on, label="capture")]


def route_schedule(J, duration, theta_A=0.0, theta_B=0.0):
    """A recombiner pulse alone, with capture off — routing calibration."""
    on = Controls(J_A=J, J_B=J, theta_A=theta_A, theta_B=theta_B)
    return [Segment(duration, on, label="route")]


def mixed_schedule(K, J, duration_per_segment=0.4, n_segments=5, seed=0):
    """A deliberately busy multi-segment schedule for the inversion test.

    Every control varies from segment to segment — capture sign, detuning,
    reference projection, recombiner strength and phase — so that the global
    inverse is exercised against a genuinely time-dependent H(t) rather than a
    single constant one.  The values come from a fixed seed and are reproducible.
    """
    rng = np.random.default_rng(seed)
    schedule = []
    for n in range(n_segments):
        controls = Controls(
            K_A=float(rng.uniform(-1.0, 1.0)) * K,
            K_B=float(rng.uniform(-1.0, 1.0)) * K,
            delta_A=float(rng.uniform(-0.5, 0.5)),
            delta_B=float(rng.uniform(-0.5, 0.5)),
            u_A=float(rng.uniform(-1.0, 1.0)),
            u_B=float(rng.uniform(-1.0, 1.0)),
            J_A=float(rng.uniform(0.0, 1.0)) * J,
            J_B=float(rng.uniform(0.0, 1.0)) * J,
            theta_A=float(rng.uniform(0.0, 2.0 * np.pi)),
            theta_B=float(rng.uniform(0.0, 2.0 * np.pi)),
        )
        schedule.append(Segment(duration_per_segment, controls,
                                label=f"mixed[{n}]"))
    return schedule
