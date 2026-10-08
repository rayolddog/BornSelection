#!/usr/bin/env python3
"""
verify.py — automated checks for the 576-state no-bath generator.

Run:

    python -m first_mark_two_absorber.verify
    python -m first_mark_two_absorber.verify --verbose
    python -m first_mark_two_absorber.verify --demo

Every check reports a measured numerical residual and the tolerance it was held
to, so a "pass" is auditable rather than a bare assertion.  The script exits 0
only if all checks pass.

What these checks establish: that the operators are the ones the plan specifies,
that the evolution is unitary, that the no-bath capture reproduces the analytic
two-level exchange, that A and B are treated identically, and that the exact
global inverse returns the complete 576-component state.  What they do NOT
establish: anything about Born weights, outcome selection, irreversibility, or
which alternative is realised.  See README.md.
"""
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

if __package__ in (None, ""):                      # allow direct script execution
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import first_mark_two_absorber                 # noqa: F401
    __package__ = "first_mark_two_absorber"

from . import controls as ctl
from . import model as mdl
from . import observables as obs
from .simulate import (final_state, phase_aligned_error, propagate,
                       segment_propagator, segment_propagator_expm)


@dataclass
class Result:
    name: str
    passed: bool
    measured: float
    tolerance: float
    detail: str = ""


CHECKS = []


def check(name):
    """Register a check function; it returns a Result."""
    def wrap(fn):
        CHECKS.append((name, fn))
        return fn
    return wrap


def _herm_residual(H):
    return float(np.max(np.abs(H - H.conj().T)))


def _busy_controls(seed=3):
    rng = np.random.default_rng(seed)
    return ctl.Controls(
        K_A=float(rng.uniform(-1, 1)), K_B=float(rng.uniform(-1, 1)),
        delta_A=float(rng.uniform(-1, 1)), delta_B=float(rng.uniform(-1, 1)),
        u_A=float(rng.uniform(-1, 1)), u_B=float(rng.uniform(-1, 1)),
        J_A=float(rng.uniform(0, 1)), J_B=float(rng.uniform(0, 1)),
        theta_A=float(rng.uniform(0, 2 * np.pi)),
        theta_B=float(rng.uniform(0, 2 * np.pi)),
    )


def _full_params():
    """Parameters with the bath couplings switched ON.

    Used only to exercise operator construction and the exact global inverse
    against a fully populated generator.  No physical conclusion about phase
    history or record formation is drawn from these runs.
    """
    return mdl.Params(chi_A=0.6, chi_B=0.45, G_A=0.5, G_B=0.5)


def _no_bath_params():
    """The no-bath setting: phase-history and mark/record couplings disabled."""
    return mdl.Params(chi_A=0.0, chi_B=0.0, G_A=0.0, G_B=0.0)


def _symmetrize(schedule):
    """Copy every arm-A control value onto arm B."""
    out = []
    for seg in schedule:
        c = seg.controls
        sym = ctl.Controls(K_A=c.K_A, K_B=c.K_A,
                           delta_A=c.delta_A, delta_B=c.delta_A,
                           u_A=c.u_A, u_B=c.u_A,
                           J_A=c.J_A, J_B=c.J_A,
                           theta_A=c.theta_A, theta_B=c.theta_A)
        out.append(ctl.Segment(seg.duration, sym, seg.sign, seg.label))
    return out


# ---------------------------------------------------------------------------
# 1. Structural checks: shapes, tensor placement, resonance relations
# ---------------------------------------------------------------------------
@check("operator shapes and tensor order")
def check_shapes():
    """Every operator is 576x576 and acts on the subsystem it names.

    Placement is verified against an independent implementation: applying the
    single-subsystem matrix along the corresponding axis of the reshaped state
    with `np.tensordot`, which never touches the Kronecker-product code path.
    """
    model = mdl.Model(_full_params())
    rng = np.random.default_rng(11)
    psi = rng.normal(size=mdl.DIM) + 1j * rng.normal(size=mdl.DIM)
    psi /= np.linalg.norm(psi)

    cases = [
        ("n_field", mdl._n_field, mdl.AXIS["field"]),
        ("S_plus_A", mdl._S_plus, mdl.AXIS["absorber_A"]),
        ("T_B", mdl._T, mdl.AXIS["absorber_B"]),
        ("Z_phase_A", mdl._Z, mdl.AXIS["phase_A"]),
        ("r_plus_B", mdl._r_plus, mdl.AXIS["record_B"]),
        ("n_r_A", mdl._n_r, mdl.AXIS["record_A"]),
    ]
    worst = 0.0
    detail = []
    for name, local, axis in cases:
        embedded = mdl.kron_op({axis: local})
        if embedded.shape != (mdl.DIM, mdl.DIM):
            return Result("operator shapes and tensor order", False, np.inf, 0.0,
                          f"{name} has shape {embedded.shape}")
        got = (embedded @ psi).reshape(mdl.DIMS)
        t = psi.reshape(mdl.DIMS)
        want = np.moveaxis(np.tensordot(local, t, axes=([1], [axis])), 0, axis)
        worst = max(worst, float(np.max(np.abs(got - want))))
        detail.append(name)

    # Two-subsystem terms must also land on the right pair of axes.
    for arm in ("A", "B"):
        ax = mdl.ARM_AXES[arm]
        a_i, _ = mdl._ARM_FIELD[arm]
        built = mdl.kron_op({mdl.AXIS["field"]: a_i, ax["absorber"]: mdl._S_plus})
        t = psi.reshape(mdl.DIMS)
        step = np.moveaxis(np.tensordot(a_i, t, axes=([1], [mdl.AXIS["field"]])),
                           0, mdl.AXIS["field"])
        step = np.moveaxis(np.tensordot(mdl._S_plus, step,
                                        axes=([1], [ax["absorber"]])),
                           0, ax["absorber"])
        worst = max(worst, float(np.max(np.abs((built @ psi).reshape(mdl.DIMS) - step))))

    shapes_ok = all(op.shape == (mdl.DIM, mdl.DIM)
                    for op in model.named_terms(_busy_controls()).values())
    if not shapes_ok:
        return Result("operator shapes and tensor order", False, np.inf, 0.0,
                      "a named term is not 576x576")
    return Result("operator shapes and tensor order", worst < 1e-12, worst, 1e-12,
                  f"dim={mdl.DIM}, order={'/'.join(mdl.SUBSYSTEMS)}")


@check("resonance relations")
def check_resonance():
    """E_e = hbar*omega and E_e = E_m + hbar*Omega_r, enforced by construction."""
    p = _full_params()
    r1 = abs(p.E_e - p.hbar * p.omega)
    r2 = abs(p.E_e - (p.E_m + p.hbar * p.Omega_r))
    worst = max(r1, r2)
    return Result("resonance relations", worst < 1e-14, worst, 1e-14,
                  f"E_e={p.E_e:.6f}, E_m={p.E_m:.6f}, hbar*Omega_r={p.hbar*p.Omega_r:.6f}")


@check("hermiticity of every generator term")
def check_hermiticity():
    """Each named term and the total H(t) are Hermitian to machine precision."""
    worst = 0.0
    where = "none (every term exactly Hermitian)"
    for params in (_no_bath_params(), _full_params()):
        model = mdl.Model(params)
        for seed in range(6):
            c = _busy_controls(seed)
            for name, op in model.named_terms(c).items():
                res = _herm_residual(op)
                if res > worst:
                    worst, where = res, name
            for sign in (1.0, -1.0):
                res = _herm_residual(model.hamiltonian(c, sign))
                if res > worst:
                    worst, where = res, f"H_total(sign={sign:+.0f})"
    return Result("hermiticity of every generator term", worst < 1e-13, worst, 1e-13,
                  f"worst term: {where}")


@check("propagator unitarity and expm cross-check")
def check_propagator():
    """exp(-iHdt) is unitary, and the eigh path agrees with scipy.linalg.expm."""
    model = mdl.Model(_full_params())
    H = model.hamiltonian(_busy_controls(2))
    U = segment_propagator(H, 0.73, model.params.hbar)
    unitarity = float(np.max(np.abs(U.conj().T @ U - np.eye(mdl.DIM))))
    U2 = segment_propagator_expm(H, 0.73, model.params.hbar)
    agreement = float(np.max(np.abs(U - U2)))
    worst = max(unitarity, agreement)
    return Result("propagator unitarity and expm cross-check", worst < 1e-11,
                  worst, 1e-11,
                  f"unitarity={unitarity:.2e}, eigh-vs-expm={agreement:.2e}")


# ---------------------------------------------------------------------------
# 2. Conservation checks
# ---------------------------------------------------------------------------
@check("norm conservation")
def check_norm():
    """The state norm stays at one through a busy multi-segment schedule."""
    worst = 0.0
    for params in (_no_bath_params(), _full_params()):
        model = mdl.Model(params)
        schedule = ctl.mixed_schedule(K=0.8, J=0.5, n_segments=6, seed=4)
        _, states = propagate(mdl.initial_state(0.7), schedule, model,
                              samples_per_segment=25)
        worst = max(worst, float(np.max(np.abs(np.linalg.norm(states, axis=1) - 1.0))))
    return Result("norm conservation", worst < 1e-12, worst, 1e-12,
                  "no-bath and full-coupling schedules, 150 sample times each")


@check("fixed-Hamiltonian energy conservation")
def check_energy():
    """With the controls held fixed, <H> does not drift."""
    worst = 0.0
    for params in (_no_bath_params(), _full_params()):
        model = mdl.Model(params)
        c = _busy_controls(5)
        schedule = [ctl.Segment(6.0, c, label="fixed")]
        _, states = propagate(mdl.initial_state(0.4), schedule, model,
                              samples_per_segment=120)
        energies = np.array([obs.expectation(s, model.hamiltonian(c)) for s in states])
        worst = max(worst, float(np.max(np.abs(energies - energies[0]))))
    return Result("fixed-Hamiltonian energy conservation", worst < 1e-11,
                  worst, 1e-11, "constant-control interval, 120 sample times")


@check("excitation number is conserved and N=2 states exist")
def check_excitation():
    """[H, N] = 0 and [H, sum_i(P_m_i - n_r_i)] = 0 — and N = 2 states are real.

    It would be wrong to say the state space forbids a second excitation.  It
    does not: `|vac, e_A, e_B>` and `|photon-A, e_A, g_B>` are ordinary vectors in
    the 576-dimensional space and both have N = 2.  This check builds them,
    confirms N = 2, and confirms they stay at N = 2 under evolution.

    What keeps the calibration runs at N = 1 is therefore dynamical, not
    structural: the prepared state starts in the N = 1 sector and H conserves N,
    so the trajectory cannot leave it.  That is a conservation law about
    excitation counting.  It is NOT evidence that only one alternative is
    realised.
    """
    model = mdl.Model(_full_params())
    worst = 0.0
    for seed in range(4):
        H = model.hamiltonian(_busy_controls(seed))
        for op in (model.N_excitation, model.mark_record_balance):
            comm = H @ op - op @ H
            worst = max(worst, float(np.max(np.abs(comm))))

    # N = 2 basis states exist in this space, and stay in their own sector.
    two_excitation = [
        mdl.basis_state(mdl.FIELD_VAC, mdl.ABS_E, mdl.ABS_E, 0, 0, 0, 0),
        mdl.basis_state(mdl.FIELD_A, mdl.ABS_E, mdl.ABS_G, 0, 0, 0, 0),
    ]
    schedule = ctl.mixed_schedule(K=0.8, J=0.5, n_segments=3, seed=12)
    sector_err = 0.0
    for psi in two_excitation:
        sector_err = max(sector_err, abs(obs.excitation_number(psi, model) - 2.0))
        _, states = propagate(psi, schedule, model, samples_per_segment=10)
        sector_err = max(sector_err, max(abs(obs.excitation_number(s, model) - 2.0)
                                         for s in states))
    # ... and the prepared calibration state stays at N = 1.
    psi1 = mdl.initial_state(0.7)
    _, states = propagate(psi1, schedule, model, samples_per_segment=10)
    sector_err = max(sector_err, max(abs(obs.excitation_number(s, model) - 1.0)
                                     for s in states))

    worst = max(worst, sector_err)
    return Result("excitation number is conserved and N=2 states exist",
                  worst < 1e-11, worst, 1e-11,
                  "[H,N]=0 and [H, mark/record balance]=0; two N=2 basis states "
                  "verified present and confined to the N=2 sector; the prepared "
                  "state stays at N=1 by conservation, not by truncation")


# ---------------------------------------------------------------------------
# 3. Analytic no-bath limits
# ---------------------------------------------------------------------------
@check("no-bath resonant capture matches analytic Rabi")
def check_rabi_resonant():
    """P(absorber i excited) = (1/2) sin^2(K t) for the equally split photon.

    On resonance the pair (photon in arm i, absorber i excited) is an exactly
    degenerate two-level system with coupling hbar*K.  The one-half is the input
    split, and the two arms do not talk to each other during capture.
    """
    model = mdl.Model(_no_bath_params())
    K = 0.9
    schedule = ctl.capture_schedule(K, duration=4.0 * np.pi / K)
    times, states = propagate(mdl.initial_state(0.7), schedule, model,
                              samples_per_segment=400)

    pe_A = np.array([obs.absorber_populations(s, "A")["excited"] for s in states])
    pe_B = np.array([obs.absorber_populations(s, "B")["excited"] for s in states])
    ph_A = np.array([obs.field_populations(s)["photon_A"] for s in states])
    analytic_e = 0.5 * np.sin(K * times) ** 2
    analytic_ph = 0.5 * np.cos(K * times) ** 2

    worst = max(float(np.max(np.abs(pe_A - analytic_e))),
                float(np.max(np.abs(pe_B - analytic_e))),
                float(np.max(np.abs(ph_A - analytic_ph))))
    return Result("no-bath resonant capture matches analytic Rabi",
                  worst < 1e-11, worst, 1e-11,
                  f"K={K}, two full Rabi cycles, 400 sample times")


@check("no-bath detuned capture matches generalized Rabi")
def check_rabi_detuned():
    """P_transfer = [K^2/(K^2+(D/2)^2)] sin^2(sqrt(K^2+(D/2)^2) t).

    A non-zero controlled detuning breaks the degeneracy, so this tests the
    detuning term as well as the capture term.
    """
    model = mdl.Model(_no_bath_params())
    K, D = 0.7, 1.1
    eff = np.sqrt(K ** 2 + (D / 2.0) ** 2)
    schedule = ctl.capture_schedule(K, duration=3.0 * np.pi / eff, delta=D)
    times, states = propagate(mdl.initial_state(0.0), schedule, model,
                              samples_per_segment=400)

    pe_A = np.array([obs.absorber_populations(s, "A")["excited"] for s in states])
    analytic = 0.5 * (K ** 2 / eff ** 2) * np.sin(eff * times) ** 2
    worst = float(np.max(np.abs(pe_A - analytic)))
    return Result("no-bath detuned capture matches generalized Rabi",
                  worst < 1e-11, worst, 1e-11,
                  f"K={K}, Delta={D}, effective rate={eff:.4f}")


@check("recombiner routes arm amplitude into the output mode")
def check_route():
    """With capture off, a J_A pulse moves photon-A into 'out' as sin^2(J t)."""
    model = mdl.Model(_no_bath_params())
    J = 0.6
    psi0 = mdl.basis_state(mdl.FIELD_A, mdl.ABS_G, mdl.ABS_G, 0, 0, 0, 0)
    schedule = [ctl.Segment(2.0 * np.pi / J,
                            ctl.Controls(J_A=J, theta_A=0.83), label="route")]
    times, states = propagate(psi0, schedule, model, samples_per_segment=300)
    p_out = np.array([obs.field_populations(s)["out"] for s in states])
    analytic = np.sin(J * times) ** 2
    worst = float(np.max(np.abs(p_out - analytic)))
    return Result("recombiner routes arm amplitude into the output mode",
                  worst < 1e-11, worst, 1e-11, f"J={J}, theta=0.83")


@check("two-arm recombiner: bright/dark interference and sqrt(2)J rate")
def check_two_arm_route():
    """Both arms drive one shared output mode, so routing is interferometric.

    With J_A = J_B = J the bright combination

        |b> = (exp(-i theta_A)|A> + exp(-i theta_B)|B>)/sqrt(2)

    couples to |out> with strength sqrt(2)*hbar*J — enhanced over the single-arm
    rate — while the orthogonal dark combination is annihilated by H_route and
    never leaves the arms.  For the prepared input the output population is

        P_out(t) = cos^2(D/2) * sin^2(sqrt(2) J t),   D = theta_in + theta_B - theta_A

    so a full transfer at t = pi/(2 sqrt(2) J) requires D = 0, and D = pi is dark.

    The declared phases are deliberately distinct and non-zero, so this check
    discriminates against a build that ignored theta or conjugated it: both of
    those give a different bright condition.  That discrimination is asserted
    numerically at the end rather than merely claimed.
    """
    model = mdl.Model(_no_bath_params())
    J, th_A, th_B = 0.6, 0.9, 0.3
    bright_time = 0.5 * np.pi / (np.sqrt(2.0) * J)

    def run(theta_in, duration, samples):
        schedule = [ctl.Segment(duration,
                                ctl.Controls(J_A=J, J_B=J,
                                             theta_A=th_A, theta_B=th_B),
                                label="route2")]
        times, states = propagate(mdl.initial_state(theta_in), schedule, model,
                                  samples_per_segment=samples)
        return times, np.array([obs.field_populations(s)["out"] for s in states])

    # (a) the full phase law, scanned over the input phase
    worst_law = 0.0
    for theta_in in np.linspace(-np.pi, np.pi, 13):
        times, p_out = run(theta_in, 3.0 * np.pi / (np.sqrt(2.0) * J), 200)
        D = theta_in + th_B - th_A
        analytic = np.cos(D / 2.0) ** 2 * np.sin(np.sqrt(2.0) * J * times) ** 2
        worst_law = max(worst_law, float(np.max(np.abs(p_out - analytic))))

    # (b) constructive input: complete transfer, and at the sqrt(2)J time
    theta_bright = th_A - th_B
    _, p_bright = run(theta_bright, bright_time, 1)
    bright_err = abs(p_bright[-1] - 1.0)

    # (c) destructive input: the dark mode never reaches the output at all
    _, p_dark = run(theta_bright + np.pi, 6.0 * np.pi / J, 300)
    dark_err = float(np.max(np.abs(p_dark)))

    # (d) the rate really is sqrt(2)J, not J: a single-arm-rate pulse undershoots
    single_arm_time = 0.5 * np.pi / J
    _, p_single = run(theta_bright, single_arm_time, 1)
    rate_gap = abs(p_single[-1] - np.sin(np.sqrt(2.0) * J * single_arm_time) ** 2)
    rate_separation = abs(p_single[-1] - 1.0)     # would be 0 if the rate were J

    # (e) sensitivity: a build that ignored theta, or that conjugated it, would
    #     have a different bright condition and so would fail (b).  The bright
    #     input here is D = 0; those two builds see D = theta_in and
    #     D = theta_in + theta_A - theta_B respectively.
    ignored = np.cos(theta_bright / 2.0) ** 2
    conjugated = np.cos((theta_bright + th_A - th_B) / 2.0) ** 2
    margin = min(abs(1.0 - ignored), abs(1.0 - conjugated))

    worst = max(worst_law, bright_err, dark_err, rate_gap)
    passed = worst < 1e-11 and rate_separation > 0.3 and margin > 0.05
    return Result("two-arm recombiner: bright/dark interference and sqrt(2)J rate",
                  passed, worst, 1e-11,
                  f"phase law over 13 input phases={worst_law:.2e}; "
                  f"bright P_out={p_bright[-1]:.12f}; dark max P_out={dark_err:.2e}; "
                  f"single-arm-rate pulse gives {p_single[-1]:.6f} (would be 1.0 if the "
                  f"bright rate were J, not sqrt(2)J); phase-convention margin={margin:.3f}")


@check("resonant return helper refuses non-resonant use")
def check_return_guard():
    """A detuned 'return' is genuinely broken, so the helper must not build one.

    Two things are checked.  First, the helper raises rather than handing back a
    schedule whose name no longer describes it.  Second — so the guard is
    protecting against something real and not just being pedantic — the schedule
    it refuses to build is constructed by hand and shown to fail: flipping the
    sign of K does not invert the detuning term, so the photon does not come back.
    """
    model = mdl.Model(_no_bath_params())
    K, t_dark, delta = 0.85, 1.9, 0.5

    raised = []
    for kwargs in ({"delta": delta}, {"u": 0.4}):
        try:
            ctl.capture_and_return_schedule(_no_bath_params(), K,
                                            t_dark=t_dark, **kwargs)
        except ctl.ResonantReturnError:
            raised.append(True)
        except Exception:                                    # wrong exception type
            raised.append(False)
        else:
            raised.append(False)
    guard_ok = all(raised)

    # The resonant helper itself must still work and still return the photon.
    psi0 = mdl.initial_state(0.7)
    good = ctl.capture_and_return_schedule(_no_bath_params(), K, t_dark=t_dark)
    good_err, _ = phase_aligned_error(psi0, final_state(psi0, good, model))
    assert all(seg.controls.delta_A == 0.0 and seg.controls.u_A == 0.0
               for seg in good), "resonant helper emitted a detuned segment"

    # Hand-built detuned flip-K schedule: the thing the guard refuses to make.
    t_pulse = ctl.full_capture_time(K)
    bad = [ctl.Segment(t_pulse, ctl.Controls(K_A=K, K_B=K,
                                             delta_A=delta, delta_B=delta)),
           ctl.Segment(t_dark, ctl.Controls(delta_A=delta, delta_B=delta)),
           ctl.Segment(t_pulse, ctl.Controls(K_A=-K, K_B=-K,
                                             delta_A=delta, delta_B=delta))]
    psi_bad = final_state(psi0, bad, model)
    bad_err, _ = phase_aligned_error(psi0, psi_bad)
    f = obs.field_populations(psi_bad)
    bad_returned = f["photon_A"] + f["photon_B"]

    passed = guard_ok and good_err < 1e-10 and bad_returned < 0.99
    return Result("resonant return helper refuses non-resonant use",
                  passed, good_err, 1e-10,
                  f"guard raises on delta and u: {guard_ok}; resonant return error="
                  f"{good_err:.2e}; hand-built detuned flip-K returns only "
                  f"{bad_returned:.10f} of the photon (state error {bad_err:.10f}) "
                  f"— which is why it is refused")


@check("no-bath return contract cannot be bypassed by omission")
def check_return_contract():
    """The mark coupling G cannot be smuggled past the helper by not mentioning it.

    G lives in `Params`, not in `Controls`, so nothing a caller passes for the
    controls reveals it.  While the physical context was optional, a caller could
    build the schedule without it, evolve under a non-zero G, and get no error
    even though the photon does not come back.  The context is now a required
    first argument, so there is no public call path that skips the check.

    Four things are checked:

    1. omitting the context is a `TypeError`, not a silently unguarded schedule;
    2. handing over a `Params` *or* a `Model` with non-zero G raises
       `ResonantReturnError` before any evolution happens;
    3. passing something that is not a Params/Model fails loudly, so the new
       leading argument cannot be satisfied by accident;
    4. the refused schedules genuinely fail, quantified over a G sweep.

    It also pins the boundary of the contract in the other direction: a non-zero
    `chi` with `G = 0` is *accepted*, because forcing `u = 0` makes the phase term
    vanish identically, and the photon does still return exactly.
    """
    K, t_dark = 0.85, 1.9
    psi0 = mdl.initial_state(0.7)
    outcomes = {}

    # 1. Omission is impossible: no context, no schedule.
    try:
        ctl.capture_and_return_schedule(K, t_dark=t_dark)
    except TypeError:
        outcomes["omitted context rejected"] = True
    except Exception:
        outcomes["omitted context rejected"] = False
    else:
        outcomes["omitted context rejected"] = False

    # 2. Non-zero G is rejected unconditionally, via Params or via Model.
    for label, context in (("Params", mdl.Params(G_A=0.5, G_B=0.5)),
                           ("Model", mdl.Model(mdl.Params(G_A=0.3, G_B=0.3))),
                           ("one-sided", mdl.Params(G_A=0.0, G_B=0.2))):
        try:
            ctl.capture_and_return_schedule(context, K, t_dark=t_dark)
        except ctl.ResonantReturnError:
            outcomes[f"G rejected via {label}"] = True
        except Exception:
            outcomes[f"G rejected via {label}"] = False
        else:
            outcomes[f"G rejected via {label}"] = False

    # 3. A wrong first argument cannot accidentally satisfy the requirement.
    for label, bogus in (("float", 0.85), ("None", None), ("dict", {"G_A": 0.0})):
        try:
            ctl.capture_and_return_schedule(bogus, K, t_dark=t_dark)
        except TypeError:
            outcomes[f"bogus context ({label}) rejected"] = True
        except Exception:
            outcomes[f"bogus context ({label}) rejected"] = False
        else:
            outcomes[f"bogus context ({label}) rejected"] = False

    # 4. What the guard is protecting against: sweep G and measure the failure.
    sweep = []
    t_pulse = ctl.full_capture_time(K)
    hand_built = [ctl.Segment(t_pulse, ctl.Controls(K_A=K, K_B=K)),
                  ctl.Segment(t_dark, ctl.Controls()),
                  ctl.Segment(t_pulse, ctl.Controls(K_A=-K, K_B=-K))]
    for G in (0.1, 0.3, 0.5):
        marked = mdl.Model(mdl.Params(G_A=G, G_B=G))
        psi_bad = final_state(psi0, hand_built, marked)
        err, _ = phase_aligned_error(psi0, psi_bad)
        f = obs.field_populations(psi_bad)
        sweep.append((G, f["photon_A"] + f["photon_B"], err))
    misuse_really_fails = all(returned < 0.99 for _, returned, _ in sweep)

    # 5. The contract boundary: chi != 0 with G = 0 is allowed, and still exact.
    chi_params = mdl.Params(chi_A=0.6, chi_B=0.45, G_A=0.0, G_B=0.0)
    chi_schedule = ctl.capture_and_return_schedule(chi_params, K, t_dark=t_dark)
    chi_err, _ = phase_aligned_error(
        psi0, final_state(psi0, chi_schedule, mdl.Model(chi_params)))
    outcomes["chi != 0 with G = 0 accepted"] = chi_err < 1e-10

    passed = all(outcomes.values()) and misuse_really_fails
    failures = [k for k, v in outcomes.items() if not v]
    sweep_text = "; ".join(f"G={G}: returns {r:.10f} (error {e:.10f})"
                           for G, r, e in sweep)
    return Result("no-bath return contract cannot be bypassed by omission",
                  passed, chi_err, 1e-10,
                  (f"all {len(outcomes)} contract probes hold; refused schedules "
                   f"do fail — {sweep_text}; chi!=0/G=0 return error {chi_err:.2e}")
                  if passed else f"FAILED probes: {failures}")


@check("mark/record exchange has the constructed frequency")
def check_mark_pair():
    """(e_A, record blank) <-> (m_A, record excited) exchanges as sin^2(G t).

    This is an operator-construction and energy-bookkeeping check on H_mark, run
    outside the no-bath setting solely to confirm the term does what its
    definition says.  It is not a record-formation or irreversibility result: a
    single record qubit simply oscillates, as the plan states.
    """
    G = 0.4
    model = mdl.Model(mdl.Params(G_A=G, G_B=G))
    psi0 = mdl.basis_state(mdl.FIELD_VAC, mdl.ABS_E, mdl.ABS_G, 0, 0, 0, 0)
    schedule = [ctl.Segment(2.0 * np.pi / G, ctl.Controls(), label="mark")]
    times, states = propagate(psi0, schedule, model, samples_per_segment=300)
    p_m = np.array([obs.absorber_populations(s, "A")["mark"] for s in states])
    p_r = np.array([obs.record_excitation(s, "A") for s in states])
    analytic = np.sin(G * times) ** 2
    worst = max(float(np.max(np.abs(p_m - analytic))),
                float(np.max(np.abs(p_r - analytic))))
    return Result("mark/record exchange has the constructed frequency",
                  worst < 1e-11, worst, 1e-11,
                  f"G={G}; mark population and record excitation track together")


# ---------------------------------------------------------------------------
# 4. Capture-and-return, symmetry, exact global inversion
# ---------------------------------------------------------------------------
@check("no-bath capture and coherent return")
def check_capture_return():
    """Capture the photon completely, hold, reverse the capture sign, get it back.

    Checks three things: the absorbers really do hold all the amplitude at the
    midpoint, the photon comes back into the arms, and the A/B relative phase
    and fringe visibility survive.
    """
    model = mdl.Model(_no_bath_params())
    K, theta_in = 0.85, 0.7
    psi0 = mdl.initial_state(theta_in)

    capture_only = ctl.capture_schedule(K, ctl.full_capture_time(K))
    mid = final_state(psi0, capture_only, model)
    captured = (obs.absorber_populations(mid, "A")["excited"]
                + obs.absorber_populations(mid, "B")["excited"])

    schedule = ctl.capture_and_return_schedule(model, K, t_dark=1.9)
    psi_end = final_state(psi0, schedule, model)
    err, overlap = phase_aligned_error(psi0, psi_end)
    returned = (obs.field_populations(psi_end)["photon_A"]
                + obs.field_populations(psi_end)["photon_B"])
    phase_err = abs(obs.arm_relative_phase(psi_end) - theta_in)
    vis = obs.interference_visibility(psi_end)

    worst = max(abs(captured - 1.0), err, abs(returned - 1.0),
                phase_err, abs(vis - 1.0))
    return Result("no-bath capture and coherent return", worst < 1e-10, worst, 1e-10,
                  f"captured={captured:.12f}, returned={returned:.12f}, "
                  f"|overlap|={overlap:.12f}, visibility={vis:.12f}")


@check("A/B symmetry under equal input and equal parameters")
def check_symmetry():
    """Equal inputs and identical parameters give identical A and B observables.

    Two statements are tested.

    (1) With theta_in = 0 the prepared state is itself invariant under exchanging
        the A and B halves, so the whole 576-component state must stay invariant
        at every sampled time.  This is run with the recombiner ON, i.e. against
        the complete generator.

    (2) With theta_in != 0 the state is no longer swap-invariant, but every A and
        B population must still match.  This is run with the recombiner OFF, and
        that is a physical requirement rather than a convenience: J_A and J_B
        both route into the *same* output mode, so once the recombiner is on the
        two arm amplitudes interfere there and the per-arm populations depend on
        the A/B relative phase by construction.  That dependence is the point of
        H_route; it is not an A/B asymmetry.  The measured size of that
        interference effect is reported in the detail line rather than hidden.
    """
    model = mdl.Model(mdl.Params(chi_A=0.5, chi_B=0.5, G_A=0.3, G_B=0.3))
    with_route = _symmetrize(ctl.mixed_schedule(K=0.8, J=0.4, n_segments=5, seed=7))
    no_route = _symmetrize(ctl.mixed_schedule(K=0.8, J=0.0, n_segments=5, seed=7))

    _, states = propagate(mdl.initial_state(0.0), with_route, model,
                          samples_per_segment=20)
    swap_err = max(float(np.max(np.abs(mdl.swap_AB(s) - s))) for s in states)

    def population_gap(schedule, theta_in):
        _, traj = propagate(mdl.initial_state(theta_in), schedule, model,
                            samples_per_segment=20)
        gap = 0.0
        for s in traj:
            pa, pb = obs.absorber_populations(s, "A"), obs.absorber_populations(s, "B")
            gap = max(gap, max(abs(pa[k] - pb[k]) for k in pa))
            gap = max(gap, abs(obs.record_excitation(s, "A")
                               - obs.record_excitation(s, "B")))
            gap = max(gap, abs(obs.phase_qubit_z(s, "A")
                               - obs.phase_qubit_z(s, "B")))
            f = obs.field_populations(s)
            gap = max(gap, abs(f["photon_A"] - f["photon_B"]))
        return gap

    pop_err = population_gap(no_route, 1.1)
    routed_gap = population_gap(with_route, 1.1)      # expected non-zero: interference
    worst = max(swap_err, pop_err)
    return Result("A/B symmetry under equal input and equal parameters",
                  worst < 1e-12, worst, 1e-12,
                  f"swap-invariance={swap_err:.2e}, population match={pop_err:.2e}; "
                  f"with the shared-output recombiner on, the phase-dependent "
                  f"arm-population difference is {routed_gap:.2e} (expected, not an "
                  f"asymmetry)")


@check("exact global inverse restores the full 576-state (no bath)")
def check_inverse_no_bath():
    """Forward, then H -> -H replayed backward, returns the initial state."""
    model = mdl.Model(_no_bath_params())
    schedule = ctl.mixed_schedule(K=0.9, J=0.6, n_segments=6, seed=8)
    psi0 = mdl.initial_state(0.7)
    psi_tau = final_state(psi0, schedule, model)
    psi_back = final_state(psi_tau, ctl.global_inverse(schedule), model)
    err, overlap = phase_aligned_error(psi0, psi_back)
    moved = float(np.linalg.norm(psi_tau - psi0))
    return Result("exact global inverse restores the full 576-state (no bath)",
                  err < 1e-11, err, 1e-11,
                  f"|overlap|={overlap:.14f}, forward displacement={moved:.6f}")


@check("exact global inverse restores the full 576-state (all couplings on)")
def check_inverse_full():
    """The same inverse with phase-history and mark/record couplings switched on.

    Included because the plan requires the global branch to reverse the complete
    generator, including the terms operational reversal may not touch.  No
    physical claim about those terms is made here.
    """
    model = mdl.Model(_full_params())
    schedule = ctl.mixed_schedule(K=1.1, J=0.7, n_segments=7, seed=9)
    psi0 = mdl.initial_state(1.3)
    psi_tau = final_state(psi0, schedule, model)
    psi_back = final_state(psi_tau, ctl.global_inverse(schedule), model)
    err, overlap = phase_aligned_error(psi0, psi_back)

    # Subsystem energies must come back too, not just the state vector.
    c0 = schedule[0].controls
    e0 = obs.energy_ledger(psi0, model, c0)
    e1 = obs.energy_ledger(psi_back, model, c0)
    ledger_err = max(abs(e0[k] - e1[k]) for k in e0)

    worst = max(err, ledger_err)
    return Result("exact global inverse restores the full 576-state (all couplings on)",
                  worst < 1e-11, worst, 1e-11,
                  f"|overlap|={overlap:.14f}, energy-ledger residual={ledger_err:.2e}")


# ---------------------------------------------------------------------------
# Demo output
# ---------------------------------------------------------------------------
def run_demo():
    """Print the no-bath capture-and-return trajectory in plain numbers."""
    model = mdl.Model(_no_bath_params())
    K = 0.85
    psi0 = mdl.initial_state(0.7)
    schedule = ctl.capture_and_return_schedule(model, K, t_dark=1.9)
    times, states = propagate(psi0, schedule, model, samples_per_segment=8)

    print("\nNo-bath capture and return (K = %.2f, theta_in = 0.70)" % K)
    print("  chi_A = chi_B = 0, G_A = G_B = 0  (phase-history and record couplings off)")
    print("  visibility is 'nan' while the photon sits in the absorbers, where the")
    print("  fringe visibility of the arm field is undefined rather than zero.")
    print()
    header = f"{'time':>7} {'photon A':>10} {'photon B':>10} {'exc A':>9} " \
             f"{'exc B':>9} {'norm':>10} {'visibility':>11}"
    print(header)
    print("-" * len(header))
    for t, s in zip(times, states):
        f = obs.field_populations(s)
        print(f"{t:7.3f} {f['photon_A']:10.6f} {f['photon_B']:10.6f} "
              f"{obs.absorber_populations(s, 'A')['excited']:9.6f} "
              f"{obs.absorber_populations(s, 'B')['excited']:9.6f} "
              f"{obs.norm(s):10.8f} {obs.interference_visibility(s):11.8f}")

    psi_end = states[-1]
    err, overlap = phase_aligned_error(psi0, psi_end)
    print(f"\nfinal |<Psi(0)|Psi(tau)>| = {overlap:.14f}")
    print(f"final state error up to global phase = {err:.3e}")
    print(f"recovered A/B relative phase = {obs.arm_relative_phase(psi_end):.14f} "
          f"(prepared 0.7)")
    print("\nThis is coherent capture and return in a closed unitary model. "
          "It is not\na statement about outcome selection, Born weights, or "
          "irreversibility.\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Verify the 576-state no-bath generator.")
    parser.add_argument("--verbose", action="store_true",
                        help="print the detail line for every check")
    parser.add_argument("--demo", action="store_true",
                        help="also print a capture-and-return trajectory")
    args = parser.parse_args(argv)

    print("first_mark_two_absorber — verification")
    print(f"  Hilbert space: {' x '.join(str(d) for d in mdl.DIMS)} = {mdl.DIM}")
    print(f"  subsystem order: {', '.join(mdl.SUBSYSTEMS)}")
    print(f"  numpy {np.__version__}\n")

    results = []
    width = max(len(name) for name, _ in CHECKS)
    for name, fn in CHECKS:
        result = fn()
        results.append(result)
        status = "PASS" if result.passed else "FAIL"
        print(f"[{status}] {result.name:<{width}}  "
              f"residual {result.measured:.3e}  (tol {result.tolerance:.0e})")
        if args.verbose or not result.passed:
            print(f"         {result.detail}")

    failed = [r for r in results if not r.passed]
    print(f"\n{len(results) - len(failed)}/{len(results)} checks passed.")

    if args.demo:
        run_demo()

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
