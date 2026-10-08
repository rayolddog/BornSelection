"""One-excitation photon/electron detector and classical threshold comparator.
All quantities dimensionless. Requires only NumPy. No files read from Downloads.
"""
from dataclasses import dataclass
import numpy as np

@dataclass
class Parameters:
    n: int = 64
    dt: float = .05
    gate: float = 12.
    g: float = 1.3
    capture: float = 1.2
    radiative: float = .10
    escape: float = .05
    thermal_sigma: float = 0.
    thermal_tau: float = .7
    random_coupling_phase: bool = False
    disorder: bool = False


def screen_state(kind='double', n=64, rng=None, count=1):
    """Conditional on transmission; Fraunhofer propagation on a periodic grid."""
    x = np.arange(n) - n//2
    left = (np.abs(x + 7) <= 1).astype(complex)
    right = (np.abs(x - 7) <= 1).astype(complex)
    phase = np.zeros(count)
    if kind == 'single':
        right *= 0
    elif kind == 'shifted':
        phase[:] = np.pi / 2
    elif kind == 'unequal':
        right *= .45
    elif kind == 'incoherent':
        if rng is None:
            raise ValueError('Random-phase ensemble requires a random generator')
        phase = rng.uniform(-np.pi, np.pi, count)
    aperture = left[None, :] + np.exp(1j*phase[:, None])*right[None, :]
    aperture /= np.linalg.norm(aperture, axis=1)[:, None]
    psi = np.fft.fftshift(np.fft.fft(aperture, axis=1, norm='ortho'), axes=1)
    return psi, phase


def intensity(kind='double', n=64):
    if kind == 'incoherent':
        # Translating identical apertures changes phase, not far-field intensity.
        return abs(screen_state('single', n)[0][0])**2
    return abs(screen_state(kind, n)[0][0])**2


def couplings(p):
    x = np.arange(p.n)
    g = np.full(p.n, p.g)
    delta = np.zeros(p.n)
    if p.disorder:
        g *= .35 + 1.15*(.5 + .5*np.sin(2*np.pi*x/p.n + .4))
        delta = 1.6*np.cos(2*np.pi*x/p.n - .8)
    return g, delta


def unitary(a, b, g, delta, t):
    """Exact exponential of H=[[0,g*],[g,delta]] at each site."""
    omega = np.sqrt(abs(g)**2 + delta**2/4)
    s = t*np.sinc(omega*t/np.pi)
    c = np.cos(omega*t)
    phase = np.exp(-.5j*delta*t)
    return (phase*((c+.5j*delta*s)*a - 1j*np.conj(g)*s*b),
            phase*(-1j*g*s*a + (c-.5j*delta*s)*b))


def damping_probabilities(a, b, p):
    ge = p.capture + p.radiative
    cap = abs(b)**2*(-np.expm1(-ge*p.dt))*p.capture/ge
    rad = abs(b)**2*(-np.expm1(-ge*p.dt))*p.radiative/ge
    esc = abs(a)**2*(-np.expm1(-p.escape*p.dt))
    return cap, rad, esc


def efficiency(p):
    """Deterministic unconditional absorption per unit input in each cell."""
    if p.thermal_sigma:
        raise ValueError('Deterministic reference requires zero thermal noise')
    g, delta = couplings(p)
    a, b = np.ones(p.n, complex), np.zeros(p.n, complex)
    capture = np.zeros(p.n)
    radiative = np.zeros(p.n)
    escape = np.zeros(p.n)
    for _ in range(round(p.gate/p.dt)):
        a, b = unitary(a, b, g, delta, p.dt/2)
        cp, rp, ep = damping_probabilities(a, b, p)
        capture += cp; radiative += rp; escape += ep
        a *= np.exp(-p.escape*p.dt/2)
        b *= np.exp(-(p.capture+p.radiative)*p.dt/2)
        a, b = unitary(a, b, g, delta, p.dt/2)
    remaining = abs(a)**2 + abs(b)**2
    return dict(capture=capture, radiative=radiative, escape=escape, remaining=remaining)


def avalanche(rng, count, target=256):
    """Phenomenological externally powered branching gain, then quenching.
    Each generation: every carrier independently adds a carrier with p=.65.
    The artificial quench caps collected charge at target electron units.
    """
    carriers = np.ones(count, dtype=int)
    latency = np.zeros(count)
    for _ in range(300):
        active = carriers < target
        if not np.any(active):
            break
        latency[active] += rng.exponential(.025, active.sum())
        carriers[active] += rng.binomial(carriers[active], .65)
        carriers = np.minimum(carriers, target)
    if np.any(carriers < target):
        raise RuntimeError('Avalanche failed to finish')
    return latency


def quantum_trials(kind, p, count=12000, seed=71231, batch=512):
    """Independent one-photon trajectories. Outcome 0=gate miss, 1=capture,
    2=radiative loss, 3=photon escape. Never samples a target position in advance.
    The jump probabilities implement standard quantum measurement postulates.
    """
    streams = np.random.SeedSequence(seed).spawn(5)
    jump_rng, phase_rng, thermal_rng, input_rng, gain_rng = [np.random.default_rng(s) for s in streams]
    outcome = np.zeros(count, int)
    site = np.full(count, -1, int)
    time = np.full(count, p.gate)
    phase_out = np.zeros(count)
    g0, delta0 = couplings(p)
    for start in range(0, count, batch):
        size = min(batch, count-start)
        a, phase = screen_state(kind, p.n, input_rng, size)
        b = np.zeros_like(a)
        phase_out[start:start+size] = phase
        g = np.broadcast_to(g0.astype(complex), a.shape).copy()
        if p.random_coupling_phase:
            g *= np.exp(1j*phase_rng.uniform(-np.pi, np.pi, a.shape))
        delta_noise = thermal_rng.normal(0, p.thermal_sigma, a.shape)
        ids = np.arange(start, start+size)
        rho = np.exp(-p.dt/p.thermal_tau)
        for step in range(round(p.gate/p.dt)):
            if not len(ids):
                break
            delta = delta0 + delta_noise
            a, b = unitary(a, b, g, delta, p.dt/2)
            cp, rp, ep = damping_probabilities(a, b, p)
            probs = np.concatenate([cp, rp, ep], axis=1)
            cdf = np.cumsum(probs, axis=1)
            u = jump_rng.random(len(ids))
            clicked = u < cdf[:, -1]
            if np.any(clicked):
                channel = np.sum(u[clicked, None] > cdf[clicked], axis=1)
                hit_ids = ids[clicked]
                outcome[hit_ids] = channel//p.n + 1
                site[hit_ids] = channel % p.n
                # Jump localized to the midpoint of this split time step.
                time[hit_ids] = (step+.5)*p.dt
            keep = ~clicked
            norm = np.sqrt(np.maximum(1-cdf[keep, -1], 1e-300))
            a = a[keep]*np.exp(-p.escape*p.dt/2)/norm[:, None]
            b = b[keep]*np.exp(-(p.capture+p.radiative)*p.dt/2)/norm[:, None]
            g = g[keep]; delta_noise = delta_noise[keep]; ids = ids[keep]
            a, b = unitary(a, b, g, delta[keep], p.dt/2)
            if p.thermal_sigma:
                delta_noise = rho*delta_noise + p.thermal_sigma*np.sqrt(1-rho*rho)*thermal_rng.normal(size=a.shape)
    latency = np.full(count, np.nan)
    latency[outcome == 1] = avalanche(gain_rng, np.count_nonzero(outcome == 1))
    return dict(outcome=outcome, site=site, time=time, latency=latency, input_phase=phase_out)


def threshold_trials(kind='double', n=64, count=12000, seed=9921,
                     law='exponential', scale=.02, noise=0., vacuum=False):
    """Finite pulse energy proxy. Independent thresholds; earliest crossing wins.
    potential_count counts crossings WITHOUT the global stop rule.
    Additive circular noise is phenomenological, not a QED vacuum model.
    """
    rng = np.random.default_rng(seed)
    psi, _ = screen_state(kind, n, rng, count)
    if vacuum:
        psi[:] = 0
    if noise:
        psi += noise/np.sqrt(2*n)*(rng.normal(size=psi.shape)+1j*rng.normal(size=psi.shape))
    drive = abs(psi)**2
    if law == 'fixed':
        threshold = np.full_like(drive, scale)
    elif law == 'uniform':
        threshold = rng.uniform(.5*scale, 1.5*scale, drive.shape)
    elif law == 'exponential':
        threshold = rng.exponential(scale, drive.shape)
    else:
        raise ValueError(law)
    fraction = np.divide(threshold, drive, out=np.full_like(drive, np.inf), where=drive>1e-28)
    best = fraction.min(axis=1)
    # Exact symmetry ties are resolved randomly, not by array index.
    ties = np.isclose(fraction, best[:, None], rtol=1e-12, atol=1e-14)
    tie_score = np.where(ties, rng.random(drive.shape), -1.)
    site = tie_score.argmax(axis=1)
    captured = best < -np.expm1(-8.)
    site[~captured] = -1
    time = np.full(count, 8.)
    time[captured] = -np.log1p(-best[captured])
    return dict(site=site, outcome=captured.astype(int), time=time,
                potential_count=(fraction < -np.expm1(-8.)).sum(axis=1))


def histogram(events, n):
    return np.bincount(events['site'][events['outcome'] == 1], minlength=n)


def goodness(counts, probabilities, seed=811, replicas=4000):
    """Multinomial Monte Carlo calibrated Pearson discrepancy, no asymptotic bins."""
    p = np.array(probabilities, float); p /= p.sum()
    total = int(np.sum(counts))
    if not total:
        return dict(n=0, tv=None, tv_null95=None, p=None)
    valid = p > 1e-14
    if np.any(np.asarray(counts)[~valid]):
        return dict(n=total, tv=float(.5*np.sum(abs(counts/total-p))), tv_null95=None, p=1/(replicas+1))
    p = p[valid]; p /= p.sum()
    obs = np.asarray(counts)[valid]
    expected = total*p
    score = np.sum((obs-expected)**2/expected)
    rng = np.random.default_rng(seed)
    null = rng.multinomial(total, p, size=replicas)
    scores = np.sum((null-expected)**2/expected, axis=1)
    tv = .5*np.sum(abs(obs/total-p))
    null_tv = .5*np.sum(abs(null/total-p), axis=1)
    return dict(n=total, tv=float(tv), tv_null95=float(np.quantile(null_tv, .95)),
                pearson=float(score), p=float((1+np.count_nonzero(scores>=score))/(replicas+1)))


def continuous_efficiency(p):
    """Independent continuous no-jump propagator with Gaussian quadrature."""
    g, delta = couplings(p)
    nodes, weights = np.polynomial.legendre.leggauss(160)
    ts = (nodes+1)*p.gate/2
    out = []
    for gj, dj in zip(g, delta):
        heff = np.array([[-.5j*p.escape, gj], [gj, dj-.5j*(p.capture+p.radiative)]])
        eig, v = np.linalg.eig(heff)
        c = np.linalg.solve(v, [1, 0])
        b = (np.exp(-1j*ts[:, None]*eig)*c) @ v[1, :]
        out.append(p.capture*p.gate/2*np.sum(weights*abs(b)**2))
    return np.array(out)
