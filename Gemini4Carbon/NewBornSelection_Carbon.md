# NewBornSelection_Carbon: Microscopic Event-by-Event Emergence of the Quantum Born Rule from Known Physics

**Date:** November 16, 2026  
**Workspace Artifacts:** `born_microscopic_sim.c`, `Makefile`, `exp1_discrete_4level.csv`, `exp2a_relaxation_ho.csv`, `exp2b_double_slit.csv`, `exp3_ehrenfest_vs_martingale.csv`

---

## Abstract

In textbook quantum mechanics, unitary state evolution and Born's probability postulate ($P = |\psi|^2$) stand as separate axioms: linear Schrödinger evolution deterministically propagates wave amplitudes, whereas measurements yield discrete individual events distributed according to $|\psi|^2$. We examine whether, and under what exact mathematical conditions, a microscopic event-by-event simulation built strictly from known physics can generate individual measurement outcomes whose accumulated statistics converge to the Born rule without ever sampling from $|\psi|^2$. First, we review the linearity no-go theorem for strictly closed pure-state Hilbert-space evolution and prove analytically why naive classical bistable symmetry breaking (Ehrenfest mean-field coupling to a thermal double-well pointer) fails, producing an error-function S-curve or step-function winner-takes-all distortion rather than $|\psi|^2$. Second, we formulate two rigorous, non-circular microscopic mechanisms grounded in established physics that overcome these obstructions: (1) **open-system continuous quantum measurement (quantum state diffusion)** driven solely by state-independent, zero-mean bath fluctuations (fair 50/50 Bernoulli coin flips) coupled to an overdamped macroscopic pointer, where norm-preserving back-action turns $|\langle k|\psi_t\rangle|^2$ into a bounded Doob martingale while preserving the exact linear Lindblad master equation in the ensemble; and (2) **configuration-space probability-current transport and stochastic quantum equilibration** (Madelung–de Broglie–Bohm hydrodynamics and Nelson–Guerra–Morato diffusion), where an exact subquantum $H$-theorem drives arbitrary non-Born initial distributions monotonically to $|\psi(x,t)|^2$ and deterministic phase-gradient advection subsequently accumulates single-particle impacts into multi-fringe interference patterns. Finally, we implement a self-contained ISO C11 simulation suite (`born_microscopic_sim.c`) and verify convergence to the Born rule across four numerical experiments via Pearson $\chi^2$ goodness-of-fit tests ($p = 0.3757$, $p = 0.1980$, and $p = 0.2478$).

---

## 1. Introduction

Since Max Born's 1926 analysis of quantum scattering [1], the relation between the quantum state vector $|\psi\rangle$ and experimental observation has rested on the Born rule: if an observable $\hat{A} = \sum_{k=1}^K a_k \hat{\Pi}_k$ is measured repeatedly on identically prepared systems in state $|\psi\rangle$, each individual run yields a single definite eigenvalue $a_k$, and the empirical relative frequencies converge to

$$
P(a_k) = \langle \psi | \hat{\Pi}_k | \psi \rangle = |\langle k | \psi \rangle|^2 \quad \text{(discrete)}, \qquad \rho(\mathbf{x}) = |\psi(\mathbf{x})|^2 \quad \text{(continuous)}.
$$

Symmetry and operational uniqueness theorems—most notably Gleason's theorem [2] for non-contextual probability measures on Hilbert spaces of dimension $d \ge 3$, Zurek's *envariance* (environment-assisted invariance) [3, 4], and operational reconstructions [5]—establish that $|\psi|^2$ is the *unique* probability measure compatible with Hilbert-space geometry and no-signaling linearity. However, uniqueness of the measure on the ensemble does not by itself specify a **microscopic event-generating dynamics**: how individual trials produce discrete detector clicks or pointer readings one by one using physical equations of motion, without presupposing an explicit random number draw weighted by $|\psi|^2$.

A genuine microscopic simulation of Born-rule selection must satisfy three strict criteria:
1. **Event-by-event definiteness:** Each individual run $r \in \{1, \dots, N\}$ must terminate in a single resolved macroscopic detector state (e.g., a pointer coordinate $M^{(r)}(T) \approx a_{k^*}$ or a localized screen impact $X^{(r)}(T)$), rather than a macroscopic superposition.
2. **Non-circularity:** No step of the algorithm may sample from $|\psi|^2$ (via inverse cumulative distribution functions, von Neumann rejection sampling, or state-biased coins). Any stochasticity representing thermal bath microstates, vacuum field shot noise, or initial particle placement must be completely independent of the quantum state amplitudes $c_k$ or $\psi(\mathbf{x}, 0)$.
3. **Physical consistency:** Ensemble averages must obey known linear quantum evolution (such as the Schrödinger continuity equation or the Lindblad master equation of open quantum systems), preventing superluminal signaling [6] and reproducing interference phenomena.

In Section 2, we analyze why closed-system linear unitary dynamics and naive semiclassical symmetry breaking fail Criteria 1 and 3, respectively. In Section 3, we derive the two exact physical mechanisms that satisfy all three criteria. Section 4 details the numerical architecture of our C implementation, and Section 5 presents statistical verification results across discrete multi-level readout, non-stationary continuous-space equilibration, Young's double-slit interference, and an Ehrenfest control test.

---

## 2. Fundamental Obstructions: Closed Unitarity and Naive Symmetry Breaking

### 2.1 The Closed-System Linearity Obstruction

Consider a tested system $S$ prepared in superposition $|\psi_S\rangle = \sum_{k=1}^K c_k |k\rangle$ ($\sum_k |c_k|^2 = 1$) coupled to a microscopic apparatus and thermal bath $D$ of $N_{\text{dof}} \gg 1$ degrees of freedom. Suppose on trial $r$ the apparatus starts in a pure microstate $|\Phi^{(r)}\rangle$ drawn from a physical ensemble $\mu(\Phi)$ (e.g., thermal phases or metastable spin configurations) completely independent of $\{c_k\}$. Under strictly linear closed-system unitary evolution $\hat{U}(t) = \exp(-i\hat{H}t/\hbar)$ on $\mathcal{H}_S \otimes \mathcal{H}_D$:

$$
|\Psi^{(r)}(t)\rangle = \hat{U}(t)\!\left(\sum_{k=1}^K c_k |k\rangle \otimes |\Phi^{(r)}\rangle\right) = \sum_{k=1}^K c_k \hat{U}(t)\big(|k\rangle \otimes |\Phi^{(r)}\rangle\big) \equiv \sum_{k=1}^K c_k |\chi_k^{(r)}(t)\rangle.
$$

Assume the apparatus is an ideal macroscopic detector on eigenstates: for each $|k\rangle$, the evolved pointer observable $\hat{M}$ satisfies $\langle\chi_k^{(r)}(t)|\hat{M}|\chi_k^{(r)}(t)\rangle \approx a_k$ with variance $O(N_{\text{dof}}^{-1})$, while environmental decoherence renders distinct branches orthogonal, $\langle\chi_j^{(r)}(t)|\hat{M}|\chi_k^{(r)}(t)\rangle \approx 0$ for $j \neq k$ [3, 7]. Evaluating the pointer expectation and variance on **any single pure-state realization $r$** yields:

$$
\langle\hat{M}\rangle^{(r)} = \langle\Psi^{(r)}(t)|\hat{M}|\Psi^{(r)}(t)\rangle \approx \sum_{k=1}^K |c_k|^2 a_k, \qquad \operatorname{Var}(\hat{M})^{(r)} \approx \sum_{k=1}^K |c_k|^2 a_k^2 - \left(\sum_{k=1}^K |c_k|^2 a_k\right)^2 > 0.
$$

Thus, as established by von Neumann [8] and formalized in general no-go theorems by Bassi and Ghirardi [9], strictly linear Hilbert-space evolution never selects a single branch on an individual pure-state run. In exactly solvable statistical-mechanical measurement models—most comprehensively the Curie–Weiss spin-boson model of Allahverdyan, Balian, and Nieuwenhuizen [7]—tracing out the phonon bath yields the reduced density operator

$$
\hat{\rho}_{S+M}(t_f) \approx \sum_{k=1}^K |c_k|^2 |k\rangle\langle k| \otimes \hat{\rho}_M^{(k)},
$$

explaining rapid off-diagonal truncation (decoherence) and ensemble registration of Born weights $|c_k|^2$, while leaving individual event selection beyond closed linear pure-state evolution.

### 2.2 Failure of Naive Semiclassical / Ehrenfest Symmetry Breaking

A natural temptation is to break the symmetry between branches by coupling a quantum state $|\psi\rangle = c_+|\uparrow\rangle + c_-|\downarrow\rangle$ to a classical bistable apparatus coordinate $X(t)$ evolving in a double-well potential $U(X) = -\frac{1}{2}a X^2 + \frac{1}{4}b X^4$ ($a, b > 0$) subject to Langevin thermal fluctuations:

$$
\gamma \frac{dX}{dt} = a X - b X^3 + g\,\langle\psi|\hat{\sigma}_z|\psi\rangle + \sqrt{2\gamma k_B T}\,\xi(t), \qquad X(0) = 0,
$$

where $\langle\xi(t)\xi(t')\rangle = \delta(t - t')$ and $\langle\psi|\hat{\sigma}_z|\psi\rangle = 2|c_+|^2 - 1$. Near the unstable symmetric hilltop $X = 0$, linearizing the escape dynamics ($\gamma \dot{X} \approx a X + h + \sqrt{2 D}\,\xi(t)$ with $h = g(2|c_+|^2 - 1)$ and $D = \gamma k_B T$) gives the exact solution

$$
X(t) = e^{a t / \gamma} \left[ \frac{h}{a}\big(1 - e^{-a t / \gamma}\big) + \frac{\sqrt{2D}}{\gamma} \int_0^t e^{-a s / \gamma}\,dW_s \right].
$$

For $t \gg \gamma / a$, the sign of $X(t)$ is governed by the Gaussian random variable $\frac{h}{a} + \mathcal{N}\!\left(0, \frac{D}{a\gamma}\right)$, yielding the selection probability

$$
P_{\text{Ehrenfest}}(+) = \frac{1}{2}\operatorname{erfc}\!\left(-\frac{g(2|c_+|^2 - 1)}{\sqrt{2 a k_B T}}\right).
$$

Rather than the linear Born law $P(+) = |c_+|^2$, naive symmetry breaking produces a Gaussian error function that degenerates in the low-noise macroscopic limit ($k_B T \ll g^2/a$) into a deterministic **winner-takes-all step function** $\Theta(|c_+|^2 - 1/2)$: a state with $|c_+|^2 = 0.51$ yields $+1$ nearly $100\%$ of the time (see also Landsman & Reuvers [10]). The missing physical ingredient is **norm-preserving quantum back-action** (or probability-flux conservation), which dynamically updates the system state as information is extracted by the bath.

---

## 3. Microscopic Physical Mechanisms Generating the Born Rule

### 3.1 Open-System Continuous Measurement and Martingale Collapse

In realistic detectors (homodyne cavity QED readout, quantum point contacts, or collisional atomic environments), the system interacts weakly and continuously with a Markovian stream of bath modes [11–15]. Under continuous quantum non-demolition (QND) monitoring of $\hat{A} = \sum_{k=1}^K a_k |k\rangle\langle k|$ ($[\hat{H}_0, \hat{A}] = 0$), the conditioned state $|\psi_t\rangle$ obeys the norm-preserving **Itô Stochastic Schrödinger Equation (SSE)** [11, 13–15]:

$$
d|\psi_t\rangle = \left[-\frac{i}{\hbar}\hat{H}_0\,dt - \frac{\kappa}{2}\big(\hat{A} - \langle\hat{A}\rangle_t\big)^2 dt + \sqrt{\kappa}\big(\hat{A} - \langle\hat{A}\rangle_t\big)\,dW_t\right]|\psi_t\rangle,
$$

where $\kappa > 0$ is the measurement coupling rate, $\langle\hat{A}\rangle_t = \langle\psi_t|\hat{A}|\psi_t\rangle$, and $dW_t$ is **unbiased standard Wiener noise** ($\mathbb{E}[dW_t] = 0$, $dW_t^2 = dt$) representing zero-mean vacuum shot noise of the bath field.

#### Exact Norm Preservation and Lindblad Consistency
Using the Itô multiplication rules $(dt)^2 = dt\,dW_t = 0$ and $(dW_t)^2 = dt$, the infinitesimal change in squared norm is

$$
d\langle\psi_t|\psi_t\rangle = -\kappa\big\langle(\hat{A} - \langle\hat{A}\rangle_t)^2\big\rangle_t\,dt + 2\sqrt{\kappa}\big\langle\hat{A} - \langle\hat{A}\rangle_t\big\rangle_t\,dW_t + \kappa\big\langle(\hat{A} - \langle\hat{A}\rangle_t)^2\big\rangle_t\,(dW_t)^2 = 0.
$$

Simultaneously, taking the expectation of $\hat{\rho}_t = \mathbb{E}[|\psi_t\rangle\langle\psi_t|]$ over the zero-mean bath fluctuations ($\mathbb{E}[dW_t] = 0$) cancels all state-dependent nonlinearities:

$$
d\hat{\rho}_t = \mathbb{E}\big[|d\psi_t\rangle\langle\psi_t| + |\psi_t\rangle\langle d\psi_t| + |d\psi_t\rangle\langle d\psi_t|\big] = -\frac{i}{\hbar}[\hat{H}_0, \hat{\rho}_t]\,dt - \frac{\kappa}{2}\big[\hat{A}, [\hat{A}, \hat{\rho}_t]\big]\,dt,
$$

which is the exact Gorini–Kossakowski–Sudarshan–Lindblad master equation of standard open-system quantum mechanics.

#### Martingale Proof of the Born Rule (Pearle–Gisin–Bauer–Bernard)
Let $p_k(t) = |\langle k | \psi_t \rangle|^2$. Applying Itô's lemma gives

$$
dp_k(t) = 2\sqrt{\kappa}\,p_k(t)\big(a_k - \langle\hat{A}\rangle_t\big)\,dW_t, \qquad \langle\hat{A}\rangle_t = \sum_{j=1}^K a_j p_j(t).
$$

Notice that **the deterministic drift $dt$ cancels identically** [15–17]:
1. **Bounded martingale:** $p_k(t) \in [0, 1]$, $\sum_{k=1}^K p_k(t) = 1$, and $\mathbb{E}[p_k(t)] = p_k(0) = |c_k|^2$ for all $t \ge 0$.
2. **Exponential localization:** Applying Itô's lemma to $\ln p_k(t)$ yields

$$
d\ln p_k(t) = -2\kappa\big(a_k - \langle\hat{A}\rangle_t\big)^2 dt + 2\sqrt{\kappa}\big(a_k - \langle\hat{A}\rangle_t\big)\,dW_t.
$$

Fixed points occur if and only if $p_k(t)(a_k - \langle\hat{A}\rangle_t) = 0$ for all $k$, which (for non-degenerate $a_k$) holds exclusively at the $K$ vertices $p_j(\infty) = \delta_{j, k^*}$. Near vertex $k^*$, all rival weights $j \neq k^*$ decay exponentially as $p_j(t) \propto \exp[-2\kappa(a_j - a_{k^*})^2 t]$.
3. **Emergence of $|\psi|^2$:** By Doob's martingale convergence theorem and bounded convergence:

$$
\mathbb{P}\big(\text{collapse to } |k^*\rangle\big) = \mathbb{E}\big[\mathbb{I}_{\{p_{k^*}(\infty) = 1\}}\big] = \mathbb{E}[p_{k^*}(\infty)] = p_{k^*}(0) = |\langle k^* | \psi(0) \rangle|^2.
$$

Coupling an overdamped physical apparatus pointer $M(t)$ via $\dot{M} = \gamma(\langle\hat{A}\rangle_t - M) + \sqrt{2\gamma D_{\text{ptr}}}\,\eta(t)$ drives $M(t) \to a_{k^*}$ on every single run, allowing events to be read directly from the macroscopic pointer position.

---

### 3.2 Configuration-Space Current Transport and Dynamical Quantum Equilibration

For continuous spatial observables, writing $\psi(\mathbf{x}, t) = R(\mathbf{x}, t) e^{i S(\mathbf{x}, t)/\hbar}$ in the time-dependent Schrödinger equation $i\hbar\partial_t\psi = \big(-\frac{\hbar^2}{2m}\nabla^2 + V\big)\psi$ yields Madelung's exact hydrodynamic continuity equation [18, 19]:

$$
\frac{\partial |\psi(\mathbf{x}, t)|^2}{\partial t} + \nabla \cdot \big(|\psi(\mathbf{x}, t)|^2 \mathbf{v}(\mathbf{x}, t)\big) = 0, \qquad \mathbf{v}(\mathbf{x}, t) = \frac{\mathbf{j}(\mathbf{x}, t)}{|\psi(\mathbf{x}, t)|^2} = \frac{\hbar}{m}\operatorname{Im}\!\left(\frac{\nabla\psi(\mathbf{x}, t)}{\psi(\mathbf{x}, t)}\right).
$$

#### Dynamical Relaxation from Arbitrary Non-Born Initial Conditions
To avoid circularity, one cannot simply initialize particle coordinates from $|\psi(\mathbf{x}, 0)|^2$. In **Nelson's stochastic mechanics** [20–22], a particle of mass $m$ undergoes conservative quantum Brownian motion with diffusion coefficient $D = \hbar/(2m)$ and forward drift $\mathbf{b}_+(\mathbf{x}, t) = \mathbf{v}(\mathbf{x}, t) + \mathbf{u}(\mathbf{x}, t)$, where $\mathbf{u} = \frac{\hbar}{m}\operatorname{Re}(\nabla\psi/\psi) = \frac{\hbar}{2m}\nabla\ln|\psi|^2$ is the osmotic velocity:

$$
d\mathbf{X}(t) = \frac{\hbar}{m}\Big[\operatorname{Im}\!\left(\frac{\nabla\psi(\mathbf{X}, t)}{\psi(\mathbf{X}, t)}\right) + \operatorname{Re}\!\left(\frac{\nabla\psi(\mathbf{X}, t)}{\psi(\mathbf{X}, t)}\right)\Big] dt + \sqrt{\frac{\hbar}{m}}\,d\mathbf{W}(t),
$$

driven purely by standard zero-mean Gaussian white noise $d\mathbf{W}(t)$. Using $\nabla^2\rho - \nabla\cdot(\rho\nabla\ln|\psi|^2) = \nabla\cdot\big(|\psi|^2 \nabla(\rho/|\psi|^2)\big)$, the Fokker–Planck equation for *any* ensemble density $\rho(\mathbf{x}, t)$ simplifies to

$$
\frac{\partial \rho}{\partial t} + \nabla\cdot(\rho\mathbf{v}) = \frac{\hbar}{2m}\nabla\cdot\left(|\psi|^2 \nabla\!\left(\frac{\rho}{|\psi|^2}\right)\right).
$$

Differentiating the Kullback–Leibler relative entropy (subquantum $H$-function) $D_{\text{KL}}(\rho \,\|\, |\psi|^2) = \int \rho \ln(\rho / |\psi|^2)\,d^3x$ and integrating by parts yields the exact **Subquantum $H$-Theorem** [21, 22]:

$$
\frac{d}{dt} D_{\text{KL}}\!\left(\rho(\cdot, t) \,\|\, |\psi(\cdot, t)|^2\right) = -\frac{\hbar}{2m}\int \rho(\mathbf{x}, t)\,\big|\nabla\ln\rho(\mathbf{x}, t) - \nabla\ln|\psi(\mathbf{x}, t)|^2\big|^2 d^3x \le 0.
$$

Consequently, even a singular non-Born initial distribution $\rho(\mathbf{x}, 0) = \delta(\mathbf{x} - \mathbf{x}_0)$ relaxes monotonically to $|\psi(\mathbf{x}, t)|^2$! (In $\ge 2$ spatial dimensions, Valentini and Westman [23, 24] proved a coarse-grained counterpart even for deterministic zero-noise de Broglie–Bohm trajectories $\dot{\mathbf{X}} = \mathbf{v}(\mathbf{X}, t)$ due to chaotic mixing around moving phase singularities.)

#### Equivariant Free-Flight Interference
Once relaxed at an aperture ($t = 0$), turning off noise ($d\mathbf{W} = 0$) leaves the purely deterministic guiding equation [19, 25]:

$$
\frac{d\mathbf{X}(t)}{dt} = \mathbf{v}(\mathbf{X}(t), t) = \frac{\hbar}{m}\operatorname{Im}\!\left(\frac{\nabla\psi(\mathbf{X}(t), t)}{\psi(\mathbf{X}(t), t)}\right).
$$

Because $\rho$ and $|\psi|^2$ satisfy the identical continuity equation $\partial_t f + \nabla\cdot(f\mathbf{v}) = 0$, equivariance guarantees $\rho(\mathbf{x}, t) = |\psi(\mathbf{x}, t)|^2$ for all $t > 0$ [25]. In a double-slit interferometer, each particle follows a continuous deterministic trajectory through *one* slit, yet the phase gradient of the superposed wave $\psi = (\psi_1 + \psi_2)/\sqrt{\mathcal{N}}$ steers trajectories away from destructive-interference nodes ($\psi = 0$) and concentrates individual impacts at constructive-interference maxima.

---

## 4. Numerical Architecture and Zero-Circularity Design

Our C11 simulation (`born_microscopic_sim.c`) implements four self-contained verification experiments using the reproducible `xoshiro256**` pseudorandom generator seeded via `SplitMix64` [26].

1. **Strict Non-Circularity:** The PRNG module exposes only three primitives:
   - `rng_fair_sign()`: extracts a single bit (`>> 63`) to return a state-independent fair 50/50 sign $\xi_n \in \{-1, +1\}$.
   - `rng_normal()`: standard Box–Muller Gaussian white noise $\mathcal{N}(0, 1)$ for Brownian increments $dW_t$.
   - `rng_uniform_open()`: flat uniform variate on $(0, 1)$ used *solely* to initialize a deliberately asymmetric non-Born distribution $X_{\text{init}} \sim \text{Uniform}[0, 2]$ on the right half-space.
2. **Exact Norm-Preserving Discrete Martingale Integrator (Experiment 1):** For step size $\Delta t$ satisfying $2\sqrt{\kappa\Delta t}\max_k|a_k - \langle\hat{A}\rangle| < 1$ and $\xi_n \in \{-1, +1\}$ ($\xi_n^2 \equiv 1$), we advance amplitudes via

$$
\psi_k^{(n+1)} = e^{-i E_k \Delta t/\hbar}\sqrt{1 + 2\sqrt{\kappa\Delta t}\big(a_k - \langle\hat{A}\rangle_n\big)\xi_n}\;\psi_k^{(n)}.
$$

   Taylor expansion in $\Delta W_n = \xi_n\sqrt{\Delta t}$ matches the exact Milstein/Itô SSE $\big[1 + \sqrt{\kappa}\Delta_k \Delta W_n - \frac{1}{2}\kappa\Delta_k^2 \Delta t + O(\Delta t^{3/2})\big]\psi_k^{(n)}$, while preserving $\sum_{k=1}^K |\psi_k^{(n+1)}|^2 = 1$ and $\mathbb{E}[|\psi_k^{(n+1)}|^2 \mid \psi^{(n)}] = |\psi_k^{(n)}|^2$ **to machine precision on every step** without ad-hoc renormalization.
3. **Second-Order Weak Stochastic Heun & $O(\Delta t^2)$ Node Taming (Experiments 2A & 2B):** Because Nelson diffusion has constant additive noise ($\sigma = \sqrt{\hbar/m} = 1$), the predictor–corrector Stochastic Heun scheme achieves weak convergence order $O(\Delta t^2)$ [27]. Isolated wavefunction nodes are regularized via quadratic taming $b_{\text{eff}} = b / \sqrt{1 + (b\,\Delta t)^2}$ [28], introducing $< 10^{-6}$ drift perturbation away from nodes while strictly bounding single-step displacements at exact zeros.
4. **Exact Quadrature & Goodness-of-Fit Statistics:** Continuous Born bin probabilities $P_{\text{Born}}(b) = \int_{x_b}^{x_{b+1}} |\psi(x, t_f)|^2\,dx$ are integrated to $< 10^{-14}$ accuracy via 8-point Gauss–Legendre quadrature, and Pearson $\chi^2$ $p$-values are computed directly in C via series and Lentz continued-fraction evaluations of the regularized incomplete gamma function $Q(\nu/2, \chi^2/2)$.

---

## 5. Simulation Results and Verification

Executing `make run` (`seed = 20261116`) runs all four experiments and exports CSV datasets.

### 5.1 Experiment 1: Discrete 4-Level Quantum System + Physical Pointer Readout

A 4-level system initialized in $|0\rangle$ evolves unitarily under an asymmetric tight-binding Hamiltonian $\hat{H}_{\text{prep}}$ for $t_{\text{prep}} = 1.15$, producing a complex superposition with irrational Born weights. Continuous QND monitoring ($\kappa = 1.0$, $\Delta t = 0.0025$, $T = 25.0$) driven strictly by fair 50/50 coin flips $\xi_n = \pm 1$ is coupled to an overdamped thermal pointer $M(t)$ ($\gamma = 2.5$, $D_{\text{ptr}} = 10^{-4}$, $M(0) = 0$). Events are classified **solely from the macroscopic pointer reading $M(T)$** ($|M(T) - a_k| < 0.25$) across $N = 200{,}000$ independent trials (`exp1_discrete_4level.csv`):

| State | Eigenvalue $a_k$ | Mean Pointer $M(T)$ | Exact Born $\|c_k\|^2$ | Empirical $P_k$ | Difference | Std. Error |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| $\|0\rangle$ | $-1.50$ | $-1.5000$ | `0.201923` | `0.202450` | `+0.000527` | `0.000898` |
| $\|1\rangle$ | $-0.50$ | $-0.5001$ | `0.469560` | `0.470875` | `+0.001315` | `0.001116` |
| $\|2\rangle$ | $+0.50$ | $+0.4999$ | `0.321770` | `0.319935` | `-0.001835` | `0.001045` |
| $\|3\rangle$ | $+1.50$ | $+1.5005$ | `0.006747` | `0.006740` | `-0.000007` | `0.000183` |

- **Resolution:** `200,000 / 200,000` resolved discrete pointer events (`0` ambiguous/unresolved runs).
- **Statistical verification:** $\chi^2 = 3.1054$ ($\text{dof} = 3$), **$p = 0.3757$**, $L_1 = 0.003683$.

---

### 5.2 Experiment 2A: Non-Equilibrium Relaxation from a Delta Point in a 3-Mode Harmonic Oscillator

We consider a non-stationary coherent superposition of the lowest three eigenstates of a 1D quantum harmonic oscillator ($\hbar = m = \omega = 1$):

$$
\psi(x, t) = \sqrt{0.45}\,\phi_0(x)e^{-i t/2} + \sqrt{0.35}\,e^{0.8 i}\phi_1(x)e^{-3 i t/2} + \sqrt{0.20}\,e^{-0.5 i}\phi_2(x)e^{-5 i t/2}.
$$

All $N = 100{,}000$ particles are initialized at the **single point** $X(0) = +1.50$ ($\rho(x, 0) = \delta(x - 1.50)$) and evolved to $t_f = 10.0$ ($\Delta t = 5 \times 10^{-4}$). Table 2 tracks the monotonic decay of the KL divergence and $L_1$ distance to the finite-sample shot-noise floor $\text{dof}/(2N) \approx 0.00016$:

| Time $t$ | $D_{\text{KL}}\big(\rho_{\text{emp}}(t) \,\|\, \|\psi(\cdot, t)\|^2\big)$ | $L_1$ Total Variation Distance |
| :---: | :---: | :---: |
| `0.10` | `0.781734` | `1.072129` |
| `0.50` | `0.289836` | `0.578226` |
| `1.50` | `0.096581` | `0.282558` |
| `3.00` | `0.010955` | `0.123347` |
| `6.00` | `0.000703` | `0.029266` |
| `10.00` | `0.000187` | `0.013393` |

At $t_f = 10.0$, $|\psi(x, 10)|^2$ exhibits an asymmetric bimodal profile (`exp2a_relaxation_ho.csv`): a primary peak at $x = -0.50$ (Born density `0.5483` vs. empirical `0.5394`), a destructive near-node trough at $x = +0.70$ (Born `0.0111` vs. empirical `0.0112`), and a secondary lobe at $x = +1.50$ (Born `0.1370` vs. empirical `0.1420`).
- **Statistical verification:** $\chi^2 = 37.4235$ ($\text{dof} = 31$), **$p = 0.1980$**.

---

### 5.3 Experiment 2B: Young's Double-Slit Single-Particle Interference Accumulation

Two Gaussian slits at $x_{1,2} = \mp d/2$ ($d = 2.20$, $\sigma_0 = 0.38$, $\hbar = m = 1$) evolve under the exact free-particle Schrödinger propagator:

$$
\psi(x, t) = \frac{(2\pi\sigma_0^2)^{-1/4}}{\sqrt{\mathcal{N}\,z(t)}}\left[\exp\!\left(-\frac{(x + d/2)^2}{4\sigma_0^2 z(t)}\right) + \exp\!\left(-\frac{(x - d/2)^2}{4\sigma_0^2 z(t)}\right)\right], \qquad z(t) = 1 + \frac{i t}{2\sigma_0^2}.
$$

1. **Asymmetric Non-Born Initialization & Equilibration:** All $N = 100{,}000$ particles start uniformly in $[0.0, 2.0]$ ($0.000\%$ left, $100.000\%$ right). Aperture equilibration over $\tau_{\text{eq}} = 48.0$ relaxes the populations across the central barrier to $49.574\%$ left / $50.426\%$ right without sampling $|\psi|^2$.
2. **Deterministic Madelung Flight:** All random noise is turned off ($dW = 0$), and each particle is advected deterministically via RK4 integration of $\dot{X} = \operatorname{Im}(\partial_x\psi / \psi)$ to $T_{\text{screen}} = 3.60$ (`exp2b_double_slit.csv`). Representative screen pixels illustrate crisp constructive and destructive interference fringes:

| Screen Coordinate $x$ | Physical Feature | Exact Born Density $\|\psi(x, T)\|^2$ | Empirical Density | Particle Count (bin $\Delta x = 0.5$) |
| :---: | :--- | :---: | :---: | :---: |
| `-8.25` | Secondary bright fringe ($m = -1$) | `0.02444` | `0.02454` | `1227` |
| `-5.25` | Destructive dark trough | `0.00164` | `0.00154` | `77` |
| `-0.25` | Central bright maximum ($m = 0$) | `0.15951` | `0.15626` | `7813` |
| `+0.25` | Central bright maximum ($m = 0$) | `0.15951` | `0.16184` | `8092` |
| `+5.25` | Destructive dark trough | `0.00164` | `0.00174` | `87` |
| `+8.25` | Secondary bright fringe ($m = +1$) | `0.02444` | `0.02436` | `1218` |

- **Statistical verification across all 48 detector bins:** $\chi^2 = 54.2715$ ($\text{dof} = 48$), **$p = 0.2478$**, $L_1 = 0.014648$.

---

### 5.4 Experiment 3: Control Comparison — Naive Ehrenfest Bistable Pointer vs. Quantum Martingale

Scanning target Born weights $|c_+|^2 \in [0.05, 0.95]$ ($40{,}000$ runs per point; `exp3_ehrenfest_vs_martingale.csv`) contrasts the naive Ehrenfest bistable pointer ($a = b = 2.0$, $g = 0.8$, $D_{\text{dw}} = 0.06$) against the 2-level norm-preserving quantum martingale coupled to a physical pointer:

| Target Born Weight $\|c_+\|^2$ | Naive Ehrenfest $P(+)$ | Ehrenfest Error | Quantum Martingale $P(+)$ | Martingale Error |
| :---: | :---: | :---: | :---: | :---: |
| `0.0500` | `0.00805` | `-0.04195` | `0.04845` | `-0.00155` |
| `0.1500` | `0.03542` | `-0.11457` | `0.14968` | `-0.00032` |
| `0.3000` | `0.16160` | `-0.13840` | `0.30238` | `+0.00238` |
| `0.5000` | `0.49925` | `-0.00075` | `0.50308` | `+0.00308` |
| `0.7000` | `0.83793` | `+0.13793` | `0.69855` | `-0.00145` |
| `0.8500` | `0.96367` | `+0.11367` | `0.84975` | `-0.00025` |
| `0.9500` | `0.99190` | `+0.04190` | `0.94978` | `-0.00022` |

Naive classical symmetry breaking exhibits systematic winner-takes-all errors up to $\pm 0.138$, whereas the unbiased quantum martingale matches $|c_+|^2$ across all points within $O(N^{-1/2})$ statistical shot noise.

---

## 6. Discussion and Conclusion

These theoretical and numerical results clarify exactly what is required for a microscopic simulation based on known physics to produce individual events obeying the Born rule:

1. **Why $|\psi|^2$ is dynamically selected:** $|\psi|^2$ is simultaneously the unique conserved density of the Schrödinger probability current ($\partial_t|\psi|^2 + \nabla\cdot\mathbf{j} = 0$) and the unique conserved martingale invariant ($\mathbb{E}[|\langle k|\psi_t\rangle|^2] = |c_k|^2$) of norm-preserving open-system diffusive monitoring.
2. **Where randomness enters without circularity:** In open-system quantum trajectories, randomness enters solely as unbiased, state-independent bath fluctuations ($\mathbb{E}[\xi_n] = 0$, literal fair coin flips); in configuration-space transport, arbitrary non-Born initial distributions relax dynamically via the subquantum $H$-theorem before deterministic wave-current steering produces interference fringes event by event.
3. **Physical interpretation:** While strictly closed, unmonitored Hilbert-space unitarity yields entanglement rather than single outcomes, coupling quantum amplitudes either to monitored open-system bath channels (operational quantum optics) or to conserved configuration-space trajectories provides rigorous, non-circular microscopic algorithms that reproduce individual quantum events and converge quantitatively to $|\psi|^2$.

---

## 7. Reproducibility

To compile and run the full verification suite from source:

```sh
make clean && make run
```

Or directly via an ISO C11 compiler (passing an optional `--seed` integer):

```sh
clang -O3 -Wall -Wextra -std=c11 born_microscopic_sim.c -lm -o born_microscopic_sim
./born_microscopic_sim --seed 20261116
```

---

## References

1. Born, M. (1926). Zur Quantenmechanik der Stoßvorgänge. *Zeitschrift für Physik*, 37(12), 863–867.
2. Gleason, A. M. (1957). Measures on the closed subspaces of a Hilbert space. *Journal of Mathematics and Mechanics*, 6(6), 885–893.
3. Zurek, W. H. (2003). Decoherence, einselection, and the quantum origins of the classical. *Reviews of Modern Physics*, 75(3), 715–775.
4. Zurek, W. H. (2005). Probabilities from entanglement, Born's rule $p_k = |\psi_k|^2$ from envariance. *Physical Review A*, 71(5), 052105.
5. Masanes, L., Galley, T. D., & Müller, M. P. (2019). The measurement postulates of quantum mechanics are operationally redundant. *Nature Communications*, 10, 1361.
6. Gisin, N. (1989). Stochastic quantum dynamics and relativity. *Helvetica Physica Acta*, 62(4), 363–371.
7. Allahverdyan, A. E., Balian, R., & Nieuwenhuizen, T. M. (2013). Understanding quantum measurement from the solution of dynamical models. *Physics Reports*, 525(1), 1–166.
8. von Neumann, J. (1932). *Mathematische Grundlagen der Quantenmechanik*. Springer, Berlin.
9. Bassi, A., & Ghirardi, G. C. (2000). A general argument against the universal validity of the superposition principle. *Physics Letters A*, 275(5–6), 373–381.
10. Landsman, N. P., & Reuvers, R. (2013). A flea on Schrödinger's cat. *Annals of Physics*, 337, 11–53.
11. Belavkin, V. P. (1989). Nondemolition measurements, nonlinear filtering and dynamic programming of quantum stochastic processes. In *Modeling and Control of Systems* (Lecture Notes in Control and Information Sciences, Vol. 121, pp. 245–265). Springer.
12. Dalibard, J., Castin, Y., & Mølmer, K. (1992). Wave-function approach to dissipative processes in quantum optics. *Physical Review Letters*, 68(5), 580–583.
13. Gisin, N., & Percival, I. C. (1992). The quantum-state diffusion model applied to open systems. *Journal of Physics A: Mathematical and General*, 25(21), 5677–5691.
14. Wiseman, H. M., & Milburn, G. J. (2009). *Quantum Measurement and Control*. Cambridge University Press.
15. Bauer, M., & Bernard, D. (2011). Convergence of repeated quantum nondemolition measurements to the Born rule. *Physical Review A*, 84(4), 044103.
16. Pearle, P. (1976). Reduction of the state vector by a nonlinear Schrödinger equation. *Physical Review D*, 13(4), 857–868.
17. Gisin, N. (1984). Quantum measurements and stochastic processes. *Physical Review Letters*, 52(19), 1657–1660.
18. Madelung, E. (1927). Quantentheorie in hydrodynamischer Form. *Zeitschrift für Physik*, 40(3–4), 322–326.
19. Bohm, D. (1952). A suggested interpretation of the quantum theory in terms of "hidden" variables. I & II. *Physical Review*, 85(2), 166–193.
20. Nelson, E. (1966). Derivation of the Schrödinger equation from Newtonian mechanics. *Physical Review*, 150(4), 1079–1085.
21. Guerra, F., & Morato, L. M. (1983). Quantization of dynamical systems and stochastic control theory. *Physical Review D*, 27(8), 1774–1786.
22. Blanchard, P., & Golin, S. (1987). Diffusion processes with singular drift fields. *Communications in Mathematical Physics*, 109(3), 421–435.
23. Valentini, A. (1991). Signal-locality, uncertainty, and the subquantum H-theorem. I. *Physics Letters A*, 156(1–2), 5–11.
24. Valentini, A., & Westman, H. (2005). Dynamical origin of quantum probabilities. *Proceedings of the Royal Society A*, 461(2053), 253–272.
25. Dürr, D., Goldstein, S., & Zanghì, N. (1992). Quantum equilibrium and the origin of absolute uncertainty. *Journal of Statistical Physics*, 67(5–6), 843–907.
26. Blackman, D., & Vigna, S. (2021). Scrambled linear pseudorandom number generators. *ACM Transactions on Mathematical Software*, 47(4), 36.
27. Kloeden, P. E., & Platen, E. (1992). *Numerical Solution of Stochastic Differential Equations*. Springer.
28. Hutzenthaler, M., Jentzen, A., & Kloeden, P. E. (2012). Strong convergence of an explicit numerical method for SDEs with nonglobally Lipschitz continuous coefficients. *Annals of Applied Probability*, 22(4), 1611–1641.
