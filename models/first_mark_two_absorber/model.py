#!/usr/bin/env python3
"""
model.py — explicit 576-dimensional operators for the two-absorber first-mark model.

Tensor order is fixed once, here, and never varies anywhere else in the package:

    field, absorber_A, absorber_B, phase_A, record_A, phase_B, record_B
      4         3           3         2         2         2         2      = 576

Every operator below is built as a full 576x576 matrix acting as the identity on
each subsystem it does not name.  Nothing in this file samples an outcome,
applies a jump, or conditions on a result: the generator is a single closed
time-dependent Hamiltonian and the evolution it produces is unitary.
"""
from dataclasses import dataclass

import numpy as np

# ---------------------------------------------------------------------------
# Subsystem ordering: the single source of truth for the whole package.
# ---------------------------------------------------------------------------
SUBSYSTEMS = ("field", "absorber_A", "absorber_B",
              "phase_A", "record_A", "phase_B", "record_B")
DIMS = (4, 3, 3, 2, 2, 2, 2)
DIM = int(np.prod(DIMS))                     # 576

AXIS = {name: k for k, name in enumerate(SUBSYSTEMS)}

# Field basis order: vacuum, photon in arm A, photon in arm B, collected output.
FIELD_VAC, FIELD_A, FIELD_B, FIELD_OUT = 0, 1, 2, 3
# Absorber basis order: ground, reversible excitation, candidate mark.
ABS_G, ABS_E, ABS_M = 0, 1, 2
# Qubit basis order: |0> / blank record, |1> / excited record.
QB_0, QB_1 = 0, 1


def _ket_bra(dim, row, col):
    """The elementary transition operator |row><col| on one subsystem."""
    op = np.zeros((dim, dim), dtype=complex)
    op[row, col] = 1.0
    return op


def kron_op(factors):
    """Embed single-subsystem operators into the full 576x576 space.

    `factors` maps a subsystem axis index to the matrix acting there; every
    unmentioned subsystem gets the identity.
    """
    out = np.ones((1, 1), dtype=complex)
    for axis, dim in enumerate(DIMS):
        out = np.kron(out, factors.get(axis, np.eye(dim, dtype=complex)))
    return out


# ---------------------------------------------------------------------------
# Basic transition operators, named exactly as in the planning artifact.
# ---------------------------------------------------------------------------
# Field: truncated single-photon transitions, not bosonic ladder operators.
_a_A = _ket_bra(4, FIELD_VAC, FIELD_A)       # remove photon from arm A
_a_B = _ket_bra(4, FIELD_VAC, FIELD_B)       # remove photon from arm B
_c_A = _ket_bra(4, FIELD_OUT, FIELD_A)       # route arm A into the output mode
_c_B = _ket_bra(4, FIELD_OUT, FIELD_B)       # route arm B into the output mode
_n_field = (_ket_bra(4, FIELD_A, FIELD_A)
            + _ket_bra(4, FIELD_B, FIELD_B)
            + _ket_bra(4, FIELD_OUT, FIELD_OUT))

# Absorber: ground -> reversible excitation -> candidate mark.
_S_plus = _ket_bra(3, ABS_E, ABS_G)          # capture into reversible excitation
_S_minus = _ket_bra(3, ABS_G, ABS_E)         # release it again
_T = _ket_bra(3, ABS_M, ABS_E)               # excitation toward candidate mark
_P_g = _ket_bra(3, ABS_G, ABS_G)
_P_e = _ket_bra(3, ABS_E, ABS_E)
_P_m = _ket_bra(3, ABS_M, ABS_M)

# Bath qubits.
_Z = np.diag([1.0, -1.0]).astype(complex)    # phase-history tag
_r_plus = _ket_bra(2, QB_1, QB_0)            # write one record excitation
_r_minus = _ket_bra(2, QB_0, QB_1)           # erase it
_n_r = _ket_bra(2, QB_1, QB_1)

# Per-arm axis lookup, so "A" and "B" code paths are literally the same code.
ARM_AXES = {
    "A": {"absorber": AXIS["absorber_A"], "phase": AXIS["phase_A"],
          "record": AXIS["record_A"]},
    "B": {"absorber": AXIS["absorber_B"], "phase": AXIS["phase_B"],
          "record": AXIS["record_B"]},
}
_ARM_FIELD = {"A": (_a_A, _c_A), "B": (_a_B, _c_B)}


@dataclass(frozen=True)
class Params:
    """Fixed physical constants of the model — not operational controls.

    The two resonance conditions of the planning artifact are enforced by
    construction rather than fitted:

        E_e = hbar * omega                (photon energy matches the excitation)
        E_e = E_m + hbar * Omega_r        (mark formation pays the record qubit)

    `chi_*` and `G_*` are the phase-history and mark/record coupling strengths.
    The no-bath setting used for the first calibration is chi = G = 0.  They live
    here, and not in `Controls`, because operational reversal is explicitly not
    permitted to flip `G_i`; only exact global inversion reverses it.
    """
    hbar: float = 1.0
    omega: float = 1.0            # photon energy is the unit of the model
    Omega_r: float = 0.35         # record-qubit gap, in units of omega
    chi_A: float = 0.0            # phase-history coupling strength at A
    chi_B: float = 0.0            # phase-history coupling strength at B
    G_A: float = 0.0              # mark/record coupling at A
    G_B: float = 0.0              # mark/record coupling at B

    @property
    def E_e(self):
        """Energy of the reversible absorber excitation."""
        return self.hbar * self.omega

    @property
    def E_m(self):
        """Energy of the candidate mark state."""
        return self.E_e - self.hbar * self.Omega_r

    def chi(self, arm):
        return self.chi_A if arm == "A" else self.chi_B

    def G(self, arm):
        return self.G_A if arm == "A" else self.G_B


class Model:
    """Builds every 576x576 operator once, then assembles H(t) on demand.

    The constant part of the generator (free field, free absorbers at zero
    detuning, free record qubits) is precomputed.  The controlled parts are kept
    as unit-coefficient Hermitian "shape" matrices so that assembling a
    Hamiltonian for a given set of control values is a handful of scalar
    multiplies rather than a rebuild.
    """

    def __init__(self, params=None):
        self.params = params if params is not None else Params()
        p = self.params

        # --- 1. free photon field: all three one-photon modes degenerate -----
        self.H_field = p.hbar * p.omega * kron_op({AXIS["field"]: _n_field})

        # --- 2/3. free absorbers (zero detuning), free bath ------------------
        self.P_e = {}
        self.P_m = {}
        self.P_g = {}
        self.n_r = {}
        self.H_absorber_free = {}
        self.H_record = {}
        for arm, ax in ARM_AXES.items():
            self.P_e[arm] = kron_op({ax["absorber"]: _P_e})
            self.P_m[arm] = kron_op({ax["absorber"]: _P_m})
            self.P_g[arm] = kron_op({ax["absorber"]: _P_g})
            self.n_r[arm] = kron_op({ax["record"]: _n_r})
            self.H_absorber_free[arm] = p.E_e * self.P_e[arm] + p.E_m * self.P_m[arm]
            self.H_record[arm] = p.hbar * p.Omega_r * self.n_r[arm]

        # The phase qubits are degenerate tags in this rung: H_phase_free = 0.
        self.H_phase_free = np.zeros((DIM, DIM), dtype=complex)
        self.H_bath = sum(self.H_record.values()) + self.H_phase_free

        # --- shape matrices for the controlled / coupled terms ---------------
        self.capture_shape = {}      # a_i S_i^+ + h.c.
        self.phase_shape = {}        # P_e_i Z_phase_i
        self.mark_shape = {}         # T_i r_i^+ + h.c.
        self.route_op = {}           # bare c_i (phase applied at assembly time)
        for arm, ax in ARM_AXES.items():
            a_i, c_i = _ARM_FIELD[arm]
            capture = kron_op({AXIS["field"]: a_i, ax["absorber"]: _S_plus})
            self.capture_shape[arm] = capture + capture.conj().T
            self.phase_shape[arm] = kron_op({ax["absorber"]: _P_e,
                                             ax["phase"]: _Z})
            mark = kron_op({ax["absorber"]: _T, ax["record"]: _r_plus})
            self.mark_shape[arm] = mark + mark.conj().T
            self.route_op[arm] = kron_op({AXIS["field"]: c_i})

        # Static generator: everything with no controllable coefficient.
        self.H_static = (self.H_field
                         + sum(self.H_absorber_free.values())
                         + self.H_bath)

        # Conserved bookkeeping operators (see README: this is excitation
        # accounting, NOT a claim about which alternative is realised).
        self.N_excitation = (kron_op({AXIS["field"]: _n_field})
                             + sum(self.P_e.values()) + sum(self.P_m.values()))
        self.mark_record_balance = (sum(self.P_m.values())
                                    - sum(self.n_r.values()))

    # -- individual named terms, for the energy ledger and for Hermiticity ----

    def H_capture(self, controls):
        """Reversible photon capture: hbar * sum_i K_i(t) [a_i S_i^+ + h.c.]."""
        p = self.params
        return p.hbar * sum(controls.K(arm) * self.capture_shape[arm]
                            for arm in ARM_AXES)

    def H_phase(self, controls):
        """Pure-dephasing phase history: hbar * sum_i (chi_i u_i(t)/2) P_e_i Z_i."""
        p = self.params
        return p.hbar * sum(0.5 * p.chi(arm) * controls.u(arm) * self.phase_shape[arm]
                            for arm in ARM_AXES)

    def H_mark(self):
        """Candidate mark plus local energy record: hbar * sum_i G_i [T_i r_i^+ + h.c.]."""
        p = self.params
        return p.hbar * sum(p.G(arm) * self.mark_shape[arm] for arm in ARM_AXES)

    def H_route(self, controls):
        """Recombiner: hbar * sum_i J_i(t) [exp(i theta_i) c_i + h.c.]."""
        p = self.params
        H = np.zeros((DIM, DIM), dtype=complex)
        for arm in ARM_AXES:
            J = controls.J(arm)
            if J == 0.0:
                continue
            phase = np.exp(1j * controls.theta(arm))
            term = phase * self.route_op[arm]
            H = H + p.hbar * J * (term + term.conj().T)
        return H

    def H_detuning(self, controls):
        """Controlled absorber detuning: hbar * sum_i Delta_i(t) P_e_i."""
        p = self.params
        return p.hbar * sum(controls.delta(arm) * self.P_e[arm] for arm in ARM_AXES)

    def named_terms(self, controls):
        """Every named piece of H(t), as a dict, for auditing and the ledger."""
        return {
            "H_field": self.H_field,
            "H_absorber_A": self.H_absorber_free["A"],
            "H_absorber_B": self.H_absorber_free["B"],
            "H_detuning": self.H_detuning(controls),
            "H_bath": self.H_bath,
            "H_capture": self.H_capture(controls),
            "H_phase": self.H_phase(controls),
            "H_mark": self.H_mark(),
            "H_route": self.H_route(controls),
        }

    def hamiltonian(self, controls, sign=1.0):
        """Assemble the complete generator H(t) for one set of control values.

        `sign = -1` gives the exact global-inversion generator -H; it is applied
        to the assembled matrix so that *every* term reverses, including the
        mark coupling and the record-qubit free energy.
        """
        H = (self.H_static
             + self.H_detuning(controls)
             + self.H_capture(controls)
             + self.H_phase(controls)
             + self.H_mark()
             + self.H_route(controls))
        return sign * H


# ---------------------------------------------------------------------------
# States
# ---------------------------------------------------------------------------
def basis_state(field, abs_A, abs_B, phase_A, record_A, phase_B, record_B):
    """A single product basis vector, indexed in the fixed subsystem order."""
    psi = np.zeros(DIMS, dtype=complex)
    psi[field, abs_A, abs_B, phase_A, record_A, phase_B, record_B] = 1.0
    return psi.reshape(-1)


def initial_state(theta_in=0.0, weight_A=None):
    """One photon coherently divided between arms A and B.

        (|photon-A> + exp(i theta_in)|photon-B>)/sqrt(2)
          (x) |g_A> (x) |g_B> (x) |+x>_A (x) |0>_A (x) |+x>_B (x) |0>_B

    Both absorbers start ready, both phase qubits in |+x> so they can acquire a
    branch-dependent phase record, both energy records blank.  `weight_A` gives
    the unequal-path variant; leave it None for the symmetric baseline.
    """
    plus_x = np.array([1.0, 1.0], dtype=complex) / np.sqrt(2.0)

    field = np.zeros(4, dtype=complex)
    if weight_A is None:
        cA = cB = 1.0 / np.sqrt(2.0)
    else:
        cA = np.sqrt(weight_A)
        cB = np.sqrt(1.0 - weight_A)
    field[FIELD_A] = cA
    field[FIELD_B] = cB * np.exp(1j * theta_in)

    ground = np.zeros(3, dtype=complex)
    ground[ABS_G] = 1.0
    blank = np.zeros(2, dtype=complex)
    blank[QB_0] = 1.0

    psi = field
    for factor in (ground, ground, plus_x, blank, plus_x, blank):
        psi = np.kron(psi, factor)
    return psi


def swap_AB(psi):
    """Exchange the A and B halves of the state: field arms, absorbers, bath."""
    t = psi.reshape(DIMS)
    t = t[[FIELD_VAC, FIELD_B, FIELD_A, FIELD_OUT], ...]     # photon-A <-> photon-B
    t = np.transpose(t, (0, 2, 1, 5, 6, 3, 4))               # absorbers and bath
    return np.ascontiguousarray(t).reshape(-1)
