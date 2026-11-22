/*
 * born_microscopic_sim.c
 *
 * Microscopic Event-by-Event Simulations Producing the Quantum Born Rule (|psi|^2)
 * Using Known Physics Dynamics Without Sampling |psi|^2
 *
 * Implements and statistically verifies three distinct experiments:
 *
 *   Experiment 1: Discrete 4-Level Quantum System + Macroscopic Apparatus Pointer
 *     - System prepared in a complex superposition via unitary Schrodinger evolution
 *       |psi(0)> = exp(-i H_prep t_prep / hbar) |0> (giving irrational |c_k|^2).
 *     - Coupled to a continuous weak QND measurement bath driven SOLELY by fair
 *       50/50 Bernoulli coin flips xi_n in {-1, +1} (zero-mean vacuum shot noise,
 *       completely independent of |psi>), plus an overdamped physical pointer M(t).
 *     - Individual events are read strictly from the macroscopic pointer M(t_final).
 *     - Verifies convergence to |c_k|^2 via Pearson chi-square goodness-of-fit.
 *
 *   Experiment 2A: Continuous-Space Non-Equilibrium Relaxation (1D Quantum Harmonic
 *                  Oscillator Superposition, Nelson Stochastic Quantum Dynamics)
 *     - Exact 3-mode time-dependent Schrodinger wavefunction psi(x, t).
 *     - ALL particles initialized at a SINGLE delta-function point X(0) = +1.50
 *       (maximally non-Born initial distribution).
 *     - Driven by Schrodinger current + osmotic drift + hbar/m Brownian diffusion.
 *     - Tracks monotonic H-theorem / KL-divergence decay and verifies final positions
 *       match the non-stationary multimodal |psi(x, t_final)|^2.
 *
 *   Experiment 2B: Young's Double-Slit Single-Particle Interference Accumulation
 *     - Particles initialized from a flat, ASYMMETRIC non-Born uniform distribution
 *       X_init ~ Uniform[0.0, 2.0] (entirely on the right half-space!).
 *     - Stage 1: Pre-flight quantum equilibration in the double-slit aperture state.
 *     - Stage 2: Purely deterministic (zero-noise) Madelung / de Broglie-Bohm
 *       probability-current advection dX/dt = Im(d_x psi / psi) during free flight
 *       to the far-field screen.
 *     - Individual screen impacts accumulate into the multi-fringe Young interference
 *       distribution |psi(x, T_screen)|^2.
 *
 *   Experiment 3: Control Comparison — Naive Ehrenfest Bistable Symmetry Breaking
 *                 ("Flea on Schrodinger's Cat") vs. Quantum Martingale Trajectory
 *     - Demonstrates quantitatively why classical double-well symmetry breaking
 *       violates the Born rule (giving an S-curve / step function), whereas the
 *       norm-preserving quantum martingale reproduces P(+) = |c_1|^2 across [0, 1].
 *
 * Compile:
 *   clang -O3 -Wall -Wextra -std=c11 born_microscopic_sim.c -lm -o born_microscopic_sim
 */

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * 1. High-Quality Reproducible PRNG: SplitMix64 + xoshiro256**
 *    (Blackman & Vigna, ACM Trans. Math. Softw. 47, 2021)
 *    Used ONLY for state-independent fair 50/50 coins, standard normal white
 *    noise, and flat non-Born initial positions. NEVER samples from |psi|^2.
 * ============================================================================ */

typedef struct {
    uint64_t s[4];
    int has_spare_normal;
    double spare_normal;
} rng_t;

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void rng_seed(rng_t *rng, uint64_t seed) {
    uint64_t sm = seed;
    for (int i = 0; i < 4; i++) {
        rng->s[i] = splitmix64(&sm);
    }
    rng->has_spare_normal = 0;
    rng->spare_normal = 0.0;
}

static inline uint64_t rotl64(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t rng_next_u64(rng_t *rng) {
    const uint64_t result = rotl64(rng->s[1] * 5ULL, 7) * 9ULL;
    const uint64_t t = rng->s[1] << 17;
    rng->s[2] ^= rng->s[0];
    rng->s[3] ^= rng->s[1];
    rng->s[1] ^= rng->s[2];
    rng->s[0] ^= rng->s[3];
    rng->s[2] ^= t;
    rng->s[3] = rotl64(rng->s[3], 45);
    return result;
}

/* Fair 50/50 symmetric sign in {-1.0, +1.0}, completely state-independent */
static inline double rng_fair_sign(rng_t *rng) {
    return (rng_next_u64(rng) >> 63) ? 1.0 : -1.0;
}

/* Uniform double in (0, 1) */
static inline double rng_uniform_open(rng_t *rng) {
    return ((rng_next_u64(rng) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/* Standard normal N(0, 1) via Box-Muller */
static inline double rng_normal(rng_t *rng) {
    if (rng->has_spare_normal) {
        rng->has_spare_normal = 0;
        return rng->spare_normal;
    }
    double u1 = rng_uniform_open(rng);
    double u2 = rng_uniform_open(rng);
    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    rng->spare_normal = r * sin(theta);
    rng->has_spare_normal = 1;
    return r * cos(theta);
}

/* ============================================================================
 * 2. Statistical Utilities: Regularized Incomplete Gamma Q(a, x) for Chi-Square
 *    p-values & 8-Point Gauss-Legendre Quadrature for Exact Bin Probabilities
 * ============================================================================ */

static double gamma_q(double a, double x) {
    if (x <= 0.0) return 1.0;
    if (x < a + 1.0) {
        /* Series representation for P(a, x), return 1 - P(a, x) */
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int n = 1; n <= 500; n++) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum) * 1e-14) break;
        }
        double p = sum * exp(-x + a * log(x) - lgamma(a));
        return 1.0 - p;
    } else {
        /* Lentz's continued fraction for Q(a, x) */
        const double eps = 1e-14;
        const double fpmin = 1e-300;
        double b = x + 1.0 - a;
        double c = 1.0 / fpmin;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= 500; i++) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (fabs(d) < fpmin) d = fpmin;
            c = b + an / c;
            if (fabs(c) < fpmin) c = fpmin;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (fabs(del - 1.0) < eps) break;
        }
        return exp(-x + a * log(x) - lgamma(a)) * h;
    }
}

static double chi2_pvalue(double chi2, int dof) {
    if (dof <= 0) return 0.0;
    return gamma_q(0.5 * (double)dof, 0.5 * chi2);
}

/* 8-point Gauss-Legendre quadrature on [a, b] */
static double quad_gl8(double (*f)(double x, const void *ctx),
                       const void *ctx, double a, double b) {
    static const double xi[4] = {
        0.183434642495649805, 0.525532409916328986,
        0.796666477413626740, 0.960289856497536232
    };
    static const double wi[4] = {
        0.362683783378361983, 0.313706645877887287,
        0.222381034453374471, 0.101228536290376259
    };
    double mid = 0.5 * (a + b);
    double half = 0.5 * (b - a);
    double sum = 0.0;
    for (int i = 0; i < 4; i++) {
        double dx = half * xi[i];
        sum += wi[i] * (f(mid - dx, ctx) + f(mid + dx, ctx));
    }
    return half * sum;
}

/* ============================================================================
 * 3. EXPERIMENT 1: Discrete 4-Level Quantum System Coupled to Continuous
 *    Unbiased Bath Fluctuations + Physical Apparatus Pointer M(t)
 * ============================================================================ */

#define K_LEVELS 4

static void matvec4(const double H[K_LEVELS][K_LEVELS],
                    const double complex psi[K_LEVELS],
                    double complex out[K_LEVELS]) {
    for (int i = 0; i < K_LEVELS; i++) {
        double complex s = 0.0;
        for (int j = 0; j < K_LEVELS; j++) {
            s += H[i][j] * psi[j];
        }
        out[i] = -I * s;
    }
}

static void prepare_4level_state(double complex psi0[K_LEVELS]) {
    /* Asymmetric 4-site Hermitian Hamiltonian (hbar = 1) */
    const double H_prep[K_LEVELS][K_LEVELS] = {
        {  0.20, -1.00,  0.00, -0.40 },
        { -1.00, -0.30, -0.80,  0.00 },
        {  0.00, -0.80,  0.50, -1.10 },
        { -0.40,  0.00, -1.10, -0.40 }
    };
    psi0[0] = 1.0 + 0.0 * I;
    psi0[1] = 0.0;
    psi0[2] = 0.0;
    psi0[3] = 0.0;

    /* Integrate i d|psi>/dt = H_prep |psi> for t_prep = 1.15 via RK4 */
    const int n_steps = 10000;
    const double dt = 1.15 / (double)n_steps;
    for (int step = 0; step < n_steps; step++) {
        double complex k1[K_LEVELS], k2[K_LEVELS], k3[K_LEVELS], k4[K_LEVELS];
        double complex tmp[K_LEVELS];

        matvec4(H_prep, psi0, k1);
        for (int i = 0; i < K_LEVELS; i++) tmp[i] = psi0[i] + 0.5 * dt * k1[i];
        matvec4(H_prep, tmp, k2);
        for (int i = 0; i < K_LEVELS; i++) tmp[i] = psi0[i] + 0.5 * dt * k2[i];
        matvec4(H_prep, tmp, k3);
        for (int i = 0; i < K_LEVELS; i++) tmp[i] = psi0[i] + dt * k3[i];
        matvec4(H_prep, tmp, k4);

        for (int i = 0; i < K_LEVELS; i++) {
            psi0[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    /* Normalize to machine precision */
    double norm2 = 0.0;
    for (int i = 0; i < K_LEVELS; i++) {
        norm2 += creal(psi0[i]) * creal(psi0[i]) + cimag(psi0[i]) * cimag(psi0[i]);
    }
    double inv_norm = 1.0 / sqrt(norm2);
    for (int i = 0; i < K_LEVELS; i++) psi0[i] *= inv_norm;
}

static void run_experiment_1(rng_t *rng) {
    printf("========================================================================\n");
    printf("EXPERIMENT 1: 4-Level Quantum System + Unbiased Continuous QND Bath\n");
    printf("              + Macroscopic Overdamped Pointer Readout M(t)\n");
    printf("========================================================================\n");

    double complex psi_init[K_LEVELS];
    prepare_4level_state(psi_init);

    double p_born[K_LEVELS];
    for (int k = 0; k < K_LEVELS; k++) {
        p_born[k] = creal(psi_init[k]) * creal(psi_init[k]) +
                    cimag(psi_init[k]) * cimag(psi_init[k]);
    }

    const double a_eig[K_LEVELS] = { -1.5, -0.5, 0.5, 1.5 };
    const double E_sys[K_LEVELS] = {  0.0,  0.3, 0.7, 1.2 };

    const int n_trials = 200000;
    const double kappa = 1.0;        /* Measurement coupling rate */
    const double gamma_ptr = 2.5;    /* Macroscopic pointer relaxation rate */
    const double D_ptr = 1.0e-4;     /* Pointer thermal diffusion */
    const double dt = 0.0025;
    const int n_steps = 10000;       /* Total time T = 25.0 */
    const double eps = 2.0 * sqrt(kappa * dt); /* = 0.10, max |eps*Delta_a| = 0.30 < 1 */
    const double ptr_noise_amp = sqrt(2.0 * gamma_ptr * D_ptr * dt);

    /* Precompute unitary phase factors per step */
    double complex phase_step[K_LEVELS];
    for (int k = 0; k < K_LEVELS; k++) {
        phase_step[k] = cos(E_sys[k] * dt) - I * sin(E_sys[k] * dt);
    }

    long counts[K_LEVELS] = {0, 0, 0, 0};
    long unresolved = 0;
    double mean_pointer_per_outcome[K_LEVELS] = {0.0, 0.0, 0.0, 0.0};

    for (int trial = 0; trial < n_trials; trial++) {
        /* Track |psi_k|^2 and complex amplitudes; note when [H_0, A] = 0,
           phase evolves unitarily while norm-preserving martingale step acts */
        double complex psi[K_LEVELS];
        for (int k = 0; k < K_LEVELS; k++) psi[k] = psi_init[k];
        double M_ptr = 0.0; /* Apparatus pointer initially at neutral 0 */

        for (int step = 0; step < n_steps; step++) {
            double p0 = creal(psi[0]) * creal(psi[0]) + cimag(psi[0]) * cimag(psi[0]);
            double p1 = creal(psi[1]) * creal(psi[1]) + cimag(psi[1]) * cimag(psi[1]);
            double p2 = creal(psi[2]) * creal(psi[2]) + cimag(psi[2]) * cimag(psi[2]);
            double p3 = creal(psi[3]) * creal(psi[3]) + cimag(psi[3]) * cimag(psi[3]);

            double mean_A = a_eig[0] * p0 + a_eig[1] * p1 +
                            a_eig[2] * p2 + a_eig[3] * p3;

            /* Fair 50/50 coin flip for vacuum fluctuation dW_t = xi * sqrt(dt) */
            double xi = rng_fair_sign(rng);

            /* Exact norm-preserving Milstein/Ito SSE step */
            psi[0] = phase_step[0] * sqrt(1.0 + eps * (a_eig[0] - mean_A) * xi) * psi[0];
            psi[1] = phase_step[1] * sqrt(1.0 + eps * (a_eig[1] - mean_A) * xi) * psi[1];
            psi[2] = phase_step[2] * sqrt(1.0 + eps * (a_eig[2] - mean_A) * xi) * psi[2];
            psi[3] = phase_step[3] * sqrt(1.0 + eps * (a_eig[3] - mean_A) * xi) * psi[3];

            /* Physical pointer relaxation + independent thermal noise */
            double eta = rng_fair_sign(rng);
            M_ptr += gamma_ptr * (mean_A - M_ptr) * dt + ptr_noise_amp * eta;
        }

        /* Read outcome STRICTLY from macroscopic pointer position M_ptr(T) */
        int detected = -1;
        for (int k = 0; k < K_LEVELS; k++) {
            if (fabs(M_ptr - a_eig[k]) < 0.25) {
                detected = k;
                break;
            }
        }
        if (detected >= 0) {
            counts[detected]++;
            mean_pointer_per_outcome[detected] += M_ptr;
        } else {
            unresolved++;
        }
    }

    int resolved_total = n_trials - (int)unresolved;
    printf("Trials: %d | Resolved macroscopic pointer events: %d | Unresolved: %ld\n\n",
           n_trials, resolved_total, unresolved);
    printf("State  Pointer a_k   Mean M(T)    Born |c_k|^2   Empirical P_k   Diff          StdErr\n");
    printf("-----  -----------   ---------    ------------   -------------   ----------    --------\n");

    double chi2 = 0.0;
    double l1_err = 0.0;
    FILE *fp = fopen("exp1_discrete_4level.csv", "w");
    if (fp) {
        fprintf(fp, "level,eigenvalue_a,mean_pointer_M,born_prob,empirical_prob,count,total_resolved\n");
    }

    for (int k = 0; k < K_LEVELS; k++) {
        double p_emp = (double)counts[k] / (double)resolved_total;
        double se = sqrt(p_born[k] * (1.0 - p_born[k]) / (double)resolved_total);
        double expected = p_born[k] * (double)resolved_total;
        double diff = p_emp - p_born[k];
        double mean_M = counts[k] > 0 ? mean_pointer_per_outcome[k] / (double)counts[k] : 0.0;

        chi2 += ((double)counts[k] - expected) * ((double)counts[k] - expected) / expected;
        l1_err += fabs(diff);

        printf(" |%d>     %+6.2f       %+7.4f      %10.6f     %10.6f     %+9.6f    %8.6f\n",
               k, a_eig[k], mean_M, p_born[k], p_emp, diff, se);
        if (fp) {
            fprintf(fp, "%d,%.4f,%.6f,%.8f,%.8f,%ld,%d\n",
                    k, a_eig[k], mean_M, p_born[k], p_emp, counts[k], resolved_total);
        }
    }
    if (fp) fclose(fp);

    int dof = K_LEVELS - 1;
    double pval = chi2_pvalue(chi2, dof);
    printf("\nGoodness-of-Fit: chi^2 = %.4f (dof = %d), p-value = %.4f, L1 distance = %.6f\n\n",
           chi2, dof, pval, l1_err);
}

/* ============================================================================
 * 4. EXPERIMENT 2A: Continuous-Space Non-Born Delta Initialization Relaxation
 *    1D Quantum Harmonic Oscillator 3-Mode Superposition (Nelson Dynamics)
 * ============================================================================ */

typedef struct {
    double complex c0;
    double complex c1;
    double complex c2;
    double t;
} ho3_ctx_t;

/* Exact |psi(x, t)|^2 for 3-mode 1D Harmonic Oscillator (hbar = m = omega = 1) */
static double ho3_born_density(double x, const void *vctx) {
    const ho3_ctx_t *ctx = (const ho3_ctx_t *)vctx;
    const double inv_pi_quarter = 0.75112554446494248286; /* pi^(-1/4) */
    const double sqrt2 = 1.41421356237309504880;
    const double inv_sqrt2 = 0.70710678118654752440;

    double complex e1 = cos(ctx->t) - I * sin(ctx->t);
    double complex e2 = cos(2.0 * ctx->t) - I * sin(2.0 * ctx->t);

    double complex Q = ctx->c0 + ctx->c1 * sqrt2 * x * e1 +
                       ctx->c2 * (2.0 * x * x - 1.0) * inv_sqrt2 * e2;
    double amp = inv_pi_quarter * exp(-0.5 * x * x);
    double q2 = creal(Q) * creal(Q) + cimag(Q) * cimag(Q);
    return amp * amp * q2;
}

/* Forward Nelson drift b_+(x, t) = v(x, t) + u(x, t) = Im(psi_x/psi) + Re(psi_x/psi) */
static inline double ho3_nelson_drift(double x, double complex c0,
                                      double complex c1_e1, double complex c2_e2,
                                      double dt) {
    const double sqrt2 = 1.41421356237309504880;
    const double inv_sqrt2 = 0.70710678118654752440;

    double complex Q = c0 + c1_e1 * (sqrt2 * x) + c2_e2 * ((2.0 * x * x - 1.0) * inv_sqrt2);
    double complex Qx = c1_e1 * sqrt2 + c2_e2 * (2.0 * sqrt2 * x);

    double denom = creal(Q) * creal(Q) + cimag(Q) * cimag(Q) + 1.0e-14;
    /* Qx / Q = (Qx * conj(Q)) / |Q|^2 */
    double re_ratio = (creal(Qx) * creal(Q) + cimag(Qx) * cimag(Q)) / denom;
    double im_ratio = (cimag(Qx) * creal(Q) - creal(Qx) * cimag(Q)) / denom;

    double b = -x + re_ratio + im_ratio;
    /* O(dt^2) quadratic taming: <1e-6 perturbation everywhere, bounded step at nodes */
    double bdt = b * dt;
    return b / sqrt(1.0 + bdt * bdt);
}

static void run_experiment_2a(rng_t *rng) {
    printf("========================================================================\n");
    printf("EXPERIMENT 2A: Continuous-Space Relaxation from Delta Point X(0) = +1.50\n");
    printf("               to Non-Stationary 3-Mode Harmonic Oscillator |psi(x,t)|^2\n");
    printf("========================================================================\n");

    const double complex c0 = sqrt(0.45);
    const double complex c1 = sqrt(0.35) * (cos(0.8) + I * sin(0.8));
    const double complex c2 = sqrt(0.20) * (cos(-0.5) + I * sin(-0.5));

    const int n_particles = 100000;
    const double x_start = 1.50; /* ALL particles start at a single point! */
    const double dt = 0.0005;
    const int n_steps = 20000;   /* Evolve from t = 0.0 to t_final = 10.0 */
    const double dt_sqrt = sqrt(dt);

    double *X = (double *)malloc((size_t)n_particles * sizeof(double));
    if (!X) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < n_particles; i++) {
        X[i] = x_start;
    }

    const int n_bins = 36;
    const double x_min = -3.6;
    const double x_max = +3.6;
    const double dx_bin = (x_max - x_min) / (double)n_bins;

    printf("Initial state: ALL %d particles placed at X(0) = +%.2f (non-Born delta!)\n",
           n_particles, x_start);
    printf("Tracking Subquantum H-function (KL divergence) decay during evolution:\n\n");
    printf("   Time t      KL Divergence D_KL(rho || |psi|^2)     L1 Total Variation\n");
    printf("   ------      ----------------------------------     ------------------\n");

    const int checkpoints[6] = { 200, 1000, 3000, 6000, 12000, 20000 };
    int cp_idx = 0;

    for (int step = 0; step < n_steps; step++) {
        double t0 = step * dt;
        double t1 = (step + 1) * dt;

        double complex c1_e1_0 = c1 * (cos(t0) - I * sin(t0));
        double complex c2_e2_0 = c2 * (cos(2.0 * t0) - I * sin(2.0 * t0));
        double complex c1_e1_1 = c1 * (cos(t1) - I * sin(t1));
        double complex c2_e2_1 = c2 * (cos(2.0 * t1) - I * sin(2.0 * t1));

        /* Stochastic Heun (2nd-order weak accuracy for additive noise) */
        for (int i = 0; i < n_particles; i++) {
            double xi = X[i];
            double dW = rng_normal(rng) * dt_sqrt;
            double b0 = ho3_nelson_drift(xi, c0, c1_e1_0, c2_e2_0, dt);
            double x_pred = xi + b0 * dt + dW;
            double b1 = ho3_nelson_drift(x_pred, c0, c1_e1_1, c2_e2_1, dt);
            X[i] = xi + 0.5 * (b0 + b1) * dt + dW;
        }

        if (cp_idx < 6 && (step + 1) == checkpoints[cp_idx]) {
            long hist[36] = {0};
            for (int i = 0; i < n_particles; i++) {
                if (X[i] >= x_min && X[i] < x_max) {
                    int b = (int)((X[i] - x_min) / dx_bin);
                    if (b >= 0 && b < n_bins) hist[b]++;
                }
            }
            ho3_ctx_t ctx = { c0, c1, c2, t1 };
            double kl = 0.0;
            double l1 = 0.0;
            for (int b = 0; b < n_bins; b++) {
                double a = x_min + b * dx_bin;
                double p_q = quad_gl8(ho3_born_density, &ctx, a, a + dx_bin);
                double p_e = (double)hist[b] / (double)n_particles;
                if (p_e > 0.0 && p_q > 1e-12) {
                    kl += p_e * log(p_e / p_q);
                }
                l1 += fabs(p_e - p_q);
            }
            printf("   %6.2f                  %10.6f                     %10.6f\n",
                   t1, kl, l1);
            cp_idx++;
        }
    }

    /* Final chi-square test at t_final = 10.0 (including combined tail bin) */
    const double t_final = n_steps * dt;
    ho3_ctx_t ctx_f = { c0, c1, c2, t_final };

    long hist[36] = {0};
    long tail_count = 0;
    for (int i = 0; i < n_particles; i++) {
        if (X[i] >= x_min && X[i] < x_max) {
            int b = (int)((X[i] - x_min) / dx_bin);
            if (b >= 0 && b < n_bins) hist[b]++;
            else tail_count++;
        } else {
            tail_count++;
        }
    }

    FILE *fp = fopen("exp2a_relaxation_ho.csv", "w");
    if (fp) {
        fprintf(fp, "bin_center,born_density,empirical_density,born_bin_prob,empirical_bin_prob,count\n");
    }

    double p_interior_sum = 0.0;
    double p_bins[36];
    for (int b = 0; b < n_bins; b++) {
        double a = x_min + b * dx_bin;
        p_bins[b] = quad_gl8(ho3_born_density, &ctx_f, a, a + dx_bin);
        p_interior_sum += p_bins[b];

        double x_mid = a + 0.5 * dx_bin;
        double rho_born = p_bins[b] / dx_bin;
        double rho_emp = ((double)hist[b] / (double)n_particles) / dx_bin;
        if (fp) {
            fprintf(fp, "%.4f,%.8f,%.8f,%.8f,%.8f,%ld\n",
                    x_mid, rho_born, rho_emp, p_bins[b],
                    (double)hist[b] / (double)n_particles, hist[b]);
        }
    }
    if (fp) fclose(fp);

    /* Merge outer sparse bins (expected < 30) into combined tail bin for exact chi^2 */
    double chi2 = 0.0;
    int active_bins = 0;
    double merged_exp = (1.0 - p_interior_sum) * (double)n_particles;
    long merged_obs = tail_count;

    for (int b = 0; b < n_bins; b++) {
        double exp_b = p_bins[b] * (double)n_particles;
        if (exp_b >= 30.0) {
            double diff = (double)hist[b] - exp_b;
            chi2 += (diff * diff) / exp_b;
            active_bins++;
        } else {
            merged_exp += exp_b;
            merged_obs += hist[b];
        }
    }
    if (merged_exp >= 5.0) {
        double diff = (double)merged_obs - merged_exp;
        chi2 += (diff * diff) / merged_exp;
        active_bins++;
    }

    int dof = active_bins - 1;
    double pval = chi2_pvalue(chi2, dof);
    printf("\nFinal Verification at t = %.2f: chi^2 = %.4f (dof = %d), p-value = %.4f\n\n",
           t_final, chi2, dof, pval);

    free(X);
}

/* ============================================================================
 * 5. EXPERIMENT 2B: Young's Double-Slit Single-Particle Interference
 *    Stage 1: Asymmetric flat Uniform[0, 2] -> slit equilibration
 *    Stage 2: Purely deterministic (zero-noise) Madelung / de Broglie-Bohm
 *             guiding equation dX/dt = Im(d_x psi / psi) to detector screen
 * ============================================================================ */

typedef struct {
    double d;        /* Slit separation */
    double sigma0;   /* Slit Gaussian width */
    double t;        /* Propagation time */
    double norm_fac; /* Normalization factor N = 2*(1 + exp(-d^2/(8*sigma0^2))) */
} dslit_ctx_t;

static double complex dslit_log_deriv(double x, const dslit_ctx_t *ctx) {
    const double s2 = ctx->sigma0 * ctx->sigma0;
    const double complex z = 1.0 + I * (ctx->t / (2.0 * s2));
    const double x1 = -0.5 * ctx->d;
    const double x2 = +0.5 * ctx->d;

    /* Compute ratio psi_2 / psi_1 stably via difference of exponents:
       -(x - x2)^2 / (4 s2 z) + (x - x1)^2 / (4 s2 z) = d * x / (2 s2 z) */
    double complex delta_arg = (ctx->d * x) / (2.0 * s2 * z);
    double re_arg = creal(delta_arg);
    double complex w1, w2;
    if (re_arg > 40.0) {
        w1 = 0.0;
        w2 = 1.0;
    } else if (re_arg < -40.0) {
        w1 = 1.0;
        w2 = 0.0;
    } else {
        double complex r = cexp(delta_arg);
        w1 = 1.0 / (1.0 + r);
        w2 = r / (1.0 + r);
    }
    return -((x - x1) * w1 + (x - x2) * w2) / (2.0 * s2 * z);
}

static double dslit_born_density(double x, const void *vctx) {
    const dslit_ctx_t *ctx = (const dslit_ctx_t *)vctx;
    const double s2 = ctx->sigma0 * ctx->sigma0;
    const double complex z = 1.0 + I * (ctx->t / (2.0 * s2));
    const double x1 = -0.5 * ctx->d;
    const double x2 = +0.5 * ctx->d;

    double pref = pow(2.0 * M_PI * s2, -0.25) / sqrt(ctx->norm_fac);
    double complex z_inv_sqrt = 1.0 / csqrt(z);

    double complex psi1 = pref * z_inv_sqrt * cexp(-((x - x1) * (x - x1)) / (4.0 * s2 * z));
    double complex psi2 = pref * z_inv_sqrt * cexp(-((x - x2) * (x - x2)) / (4.0 * s2 * z));
    double complex psi = psi1 + psi2;
    return creal(psi) * creal(psi) + cimag(psi) * cimag(psi);
}

static inline double dslit_velocity(double x, double t, double d, double sigma0) {
    dslit_ctx_t ctx = { d, sigma0, t, 1.0 };
    return cimag(dslit_log_deriv(x, &ctx));
}

static void run_experiment_2b(rng_t *rng) {
    printf("========================================================================\n");
    printf("EXPERIMENT 2B: Young's Double-Slit Interference Event-by-Event\n");
    printf("               (Asymmetric Flat Init -> Deterministic Madelung Flight)\n");
    printf("========================================================================\n");

    const double d = 2.20;
    const double sigma0 = 0.38;
    const double norm_fac = 2.0 * (1.0 + exp(-(d * d) / (8.0 * sigma0 * sigma0)));
    const double T_screen = 3.60;
    const int n_particles = 100000;

    double *X = (double *)malloc((size_t)n_particles * sizeof(double));
    if (!X) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    /* Non-Born Initialization: ALL particles placed uniformly on RIGHT side [0.0, 2.0] */
    for (int i = 0; i < n_particles; i++) {
        X[i] = 2.0 * rng_uniform_open(rng);
    }

    /* Stage 1: Aperture quantum equilibration at t = 0.
       At t = 0, d_x ln psi(x,0) = -(x - (d/2)*tanh(d*x / (4*sigma0^2))) / (2*sigma0^2).
       Integrate over tau_eq = 48.0 (>12 Kramers inter-slit barrier crossing times)
       so the 100%-right initial cloud equilibrates across both slits without sampling |psi|^2. */
    const double dtau = 0.002;
    const int n_eq_steps = 24000;
    const double dtau_sqrt = sqrt(dtau);
    const double inv_2s2 = 1.0 / (2.0 * sigma0 * sigma0);
    const double tanh_scale = d / (4.0 * sigma0 * sigma0);
    const double half_d = 0.5 * d;

    for (int step = 0; step < n_eq_steps; step++) {
        for (int i = 0; i < n_particles; i++) {
            double xi = X[i];
            double dW = rng_normal(rng) * dtau_sqrt;
            double b0 = -(xi - half_d * tanh(tanh_scale * xi)) * inv_2s2;
            double x_pred = xi + b0 * dtau + dW;
            double b1 = -(x_pred - half_d * tanh(tanh_scale * x_pred)) * inv_2s2;
            X[i] = xi + 0.5 * (b0 + b1) * dtau + dW;
        }
    }

    long left_slit = 0, right_slit = 0;
    for (int i = 0; i < n_particles; i++) {
        if (X[i] < 0.0) left_slit++;
        else right_slit++;
    }
    printf("Stage 1 Equilibration (from 0.00%% left / 100.00%% right non-Born init):\n");
    printf("  Left slit (x < 0): %.3f%% | Right slit (x >= 0): %.3f%%\n\n",
           100.0 * (double)left_slit / (double)n_particles,
           100.0 * (double)right_slit / (double)n_particles);

    /* Stage 2: 100% DETERMINISTIC free-space Madelung / de Broglie-Bohm flight
       dX/dt = j(X, t) / |psi(X, t)|^2 = Im(d_x psi / psi) via RK4 (ZERO noise!) */
    const int n_flight_steps = 1800;
    const double dt = T_screen / (double)n_flight_steps;

    for (int step = 0; step < n_flight_steps; step++) {
        double t0 = step * dt;
        double th = t0 + 0.5 * dt;
        double t1 = t0 + dt;
        for (int i = 0; i < n_particles; i++) {
            double x0 = X[i];
            double k1 = dslit_velocity(x0, t0, d, sigma0);
            double k2 = dslit_velocity(x0 + 0.5 * dt * k1, th, d, sigma0);
            double k3 = dslit_velocity(x0 + 0.5 * dt * k2, th, d, sigma0);
            double k4 = dslit_velocity(x0 + dt * k3, t1, d, sigma0);
            X[i] = x0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
    }

    /* Screen detector array: 48 pixels across [-12.0, +12.0] */
    const int n_bins = 48;
    const double x_min = -12.0;
    const double x_max = +12.0;
    const double dx_bin = (x_max - x_min) / (double)n_bins;

    long hist[48] = {0};
    long tail_count = 0;
    for (int i = 0; i < n_particles; i++) {
        if (X[i] >= x_min && X[i] < x_max) {
            int b = (int)((X[i] - x_min) / dx_bin);
            if (b >= 0 && b < n_bins) hist[b]++;
            else tail_count++;
        } else {
            tail_count++;
        }
    }

    dslit_ctx_t ctx_screen = { d, sigma0, T_screen, norm_fac };
    double p_bins[48];
    double p_interior_sum = 0.0;
    double l1_err = 0.0;

    FILE *fp = fopen("exp2b_double_slit.csv", "w");
    if (fp) {
        fprintf(fp, "x_screen,born_density,empirical_density,born_prob,empirical_prob,count\n");
    }

    for (int b = 0; b < n_bins; b++) {
        double a = x_min + b * dx_bin;
        p_bins[b] = quad_gl8(dslit_born_density, &ctx_screen, a, a + dx_bin);
        p_interior_sum += p_bins[b];
        double p_emp = (double)hist[b] / (double)n_particles;
        l1_err += fabs(p_emp - p_bins[b]);

        if (fp) {
            fprintf(fp, "%.4f,%.8f,%.8f,%.8f,%.8f,%ld\n",
                    a + 0.5 * dx_bin, p_bins[b] / dx_bin, p_emp / dx_bin,
                    p_bins[b], p_emp, hist[b]);
        }
    }
    if (fp) fclose(fp);

    /* Display ASCII profile of central interference fringes */
    printf("Screen Interference Profile at T = %.2f (central pixels, '*' = Empirical, '|' = Born):\n\n",
           T_screen);
    for (int b = 8; b < n_bins - 8; b++) {
        double x_mid = x_min + (b + 0.5) * dx_bin;
        double rho_b = p_bins[b] / dx_bin;
        double rho_e = ((double)hist[b] / (double)n_particles) / dx_bin;
        int col_b = (int)round(rho_b * 220.0);
        int col_e = (int)round(rho_e * 220.0);
        if (col_b > 48) col_b = 48;
        if (col_e > 48) col_e = 48;

        char line[52];
        memset(line, ' ', 50);
        line[50] = '\0';
        for (int c = 0; c < col_e; c++) line[c] = '*';
        if (col_b >= 0 && col_b < 50) line[col_b] = '|';
        printf(" x=%+5.2f [%6.4f vs %6.4f] |%s\n", x_mid, rho_e, rho_b, line);
    }

    /* Exact Pearson chi-square test (merging bins with expected < 25) */
    double chi2 = 0.0;
    int active_bins = 0;
    double merged_exp = (1.0 - p_interior_sum) * (double)n_particles;
    long merged_obs = tail_count;

    for (int b = 0; b < n_bins; b++) {
        double exp_b = p_bins[b] * (double)n_particles;
        if (exp_b >= 25.0) {
            double diff = (double)hist[b] - exp_b;
            chi2 += (diff * diff) / exp_b;
            active_bins++;
        } else {
            merged_exp += exp_b;
            merged_obs += hist[b];
        }
    }
    if (merged_exp >= 5.0) {
        double diff = (double)merged_obs - merged_exp;
        chi2 += (diff * diff) / merged_exp;
        active_bins++;
    }

    int dof = active_bins - 1;
    double pval = chi2_pvalue(chi2, dof);
    printf("\nDouble-Slit Goodness-of-Fit: chi^2 = %.4f (dof = %d), p-value = %.4f, L1 = %.6f\n\n",
           chi2, dof, pval, l1_err);

    free(X);
}

/* ============================================================================
 * 6. EXPERIMENT 3: Control Comparison — Why Naive Ehrenfest Double-Well
 *    Symmetry Breaking Violates Born's Rule While Quantum Martingale Succeeds
 * ============================================================================ */

static void run_experiment_3(rng_t *rng) {
    printf("========================================================================\n");
    printf("EXPERIMENT 3: Control Test across Born Weights p = |c_+|^2 in [0.05, 0.95]\n");
    printf("              Naive Ehrenfest Bistable Pointer vs. Quantum Martingale\n");
    printf("========================================================================\n\n");

    const double targets[7] = { 0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95 };
    const int n_runs = 40000;

    /* Ehrenfest double-well dX = (a X - b X^3 + g*(2p - 1))*dt + sqrt(2 D)*dW */
    const double a_dw = 2.0, b_dw = 2.0, g_dw = 0.8, D_dw = 0.06;
    const double dt_dw = 0.005;
    const int steps_dw = 800;
    const double noise_dw = sqrt(2.0 * D_dw * dt_dw);

    /* Quantum martingale 2-level dp = 2*sqrt(kappa)*p*(1-p)*dW */
    const double kappa = 1.0, dt_qm = 0.004;
    const int steps_qm = 3500;
    const double eps_qm = 2.0 * sqrt(kappa * dt_qm);

    printf("Target |c_+|^2   Naive Ehrenfest P(+)   Ehrenfest Err   Martingale P(+)   Martingale Err\n");
    printf("--------------   --------------------   -------------   ---------------   --------------\n");

    FILE *fp = fopen("exp3_ehrenfest_vs_martingale.csv", "w");
    if (fp) {
        fprintf(fp, "born_target,ehrenfest_prob,ehrenfest_error,martingale_prob,martingale_error\n");
    }

    for (int idx = 0; idx < 7; idx++) {
        double p0 = targets[idx];
        double bias = g_dw * (2.0 * p0 - 1.0);

        long wins_ehrenfest = 0;
        long wins_martingale = 0;

        for (int r = 0; r < n_runs; r++) {
            /* (a) Naive Ehrenfest double-well starting at hilltop X(0) = 0 */
            double X = 0.0;
            for (int s = 0; s < steps_dw; s++) {
                X += (a_dw * X - b_dw * X * X * X + bias) * dt_dw +
                     noise_dw * rng_fair_sign(rng);
            }
            if (X > 0.0) wins_ehrenfest++;

            /* (b) Unbiased 2-level quantum martingale trajectory + pointer */
            double p = p0;
            double M = 0.0;
            for (int s = 0; s < steps_qm; s++) {
                double mean_sz = 2.0 * p - 1.0; /* eigenvalues +1, -1 */
                double xi = rng_fair_sign(rng);
                /* Exact norm-preserving update p_{n+1} = p_n * (1 + eps*(1 - mean_sz)*xi) */
                p = p * (1.0 + eps_qm * (1.0 - mean_sz) * xi);
                M += 2.5 * (mean_sz - M) * dt_qm;
            }
            if (M > 0.0) wins_martingale++;
        }

        double p_ehr = (double)wins_ehrenfest / (double)n_runs;
        double p_mar = (double)wins_martingale / (double)n_runs;

        printf("    %6.4f              %8.5f          %+8.5f           %8.5f         %+8.5f\n",
               p0, p_ehr, p_ehr - p0, p_mar, p_mar - p0);
        if (fp) {
            fprintf(fp, "%.4f,%.6f,%.6f,%.6f,%.6f\n",
                    p0, p_ehr, p_ehr - p0, p_mar, p_mar - p0);
        }
    }
    if (fp) fclose(fp);
    printf("\n");
}

int main(int argc, char **argv) {
    uint64_t seed = 20261116ULL;
    if (argc > 2 && strcmp(argv[1], "--seed") == 0) {
        seed = (uint64_t)strtoull(argv[2], NULL, 10);
    }

    rng_t rng;
    rng_seed(&rng, seed);

    printf("Microscopic Quantum Event Simulation Suite (PRNG seed = %llu)\n\n",
           (unsigned long long)seed);

    run_experiment_1(&rng);
    run_experiment_2a(&rng);
    run_experiment_2b(&rng);
    run_experiment_3(&rng);

    printf("All CSV verification outputs written to current workspace directory.\n");
    return 0;
}
