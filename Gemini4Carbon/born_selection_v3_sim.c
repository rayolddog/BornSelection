/*
 * born_selection_v3_sim.c
 *
 * Numerical verification comparing two single-event simulation protocols
 * for a Schrödinger wave packet illuminating a row of solid-state detectors:
 *
 *   Protocol A (Copenhagen / Born Rule + RNG):
 *     1. Propagate psi_inc(x, T) via the time-dependent Schrödinger equation.
 *     2. Apply the Born rule P_m = int_m |psi_inc(x, T)|^2 dx to form a
 *        probability distribution.
 *     3. Draw a single detector outcome via inverse-CDF random sampling.
 *
 *   Protocol B (Proposed Real-Wave Surface-Interference, First-to-Tunnel,
 *               and Random-Phase Solid-State Cascade Mechanism):
 *     1. Propagate the same real wavepsi_inc(x, T) via the Schrödinger equation.
 *     2. At each near-surface atomic site (m, a) across M = 60 detectors
 *        (N_surf = 128 surface atoms per detector), compute the complex three-wave
 *        superposition of the incident wave packet, bound surface-electron wave
 *        packets with random phases, and zero-point quantum-vacuum fluctuations.
 *        Subtracting the unilluminated equilibrium bound+vacuum background gives
 *        local deposited energy with strictly unbiased interference fluctuations.
 *     3. Each surface electron has a state-independent pre-arrival sub-barrier
 *        action deficit Delta_S_{m,a} ~ Uniform(0, 1) below the orbital tunneling
 *        barrier (no sampling of |psi|^2 anywhere).
 *     4. The FIRST atomic site (m*, a*) to tunnel absorbs the entire photon
 *        quantum via wave-packet nonlocality (suppressing all other sites, g^(2)(0)=0).
 *     5. The tunneled electron triggers an inelastic cascade across solid-state
 *        bound electrons with randomly distributed phases, scrambling incident
 *        wave phase information irreversibly.
 *
 * Compile & run:
 *   clang -O3 -Wall -Wextra -std=c11 born_selection_v3_sim.c -lm -o born_selection_v3_sim
 *   ./born_selection_v3_sim
 */

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    uint64_t s[4];
} rng_t;

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void rng_seed(rng_t *rng, uint64_t seed) {
    uint64_t sm = seed;
    for (int i = 0; i < 4; i++) rng->s[i] = splitmix64(&sm);
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

static inline double rng_uniform_open(rng_t *rng) {
    return ((double)(rng_next_u64(rng) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/* Regularized upper incomplete gamma Q(a, x) for Pearson chi-square p-values */
static double gamma_q(double a, double x) {
    if (x <= 0.0) return 1.0;
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int n = 1; n <= 500; n++) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum) * 1e-14) break;
        }
        return 1.0 - sum * exp(-x + a * log(x) - lgamma(a));
    } else {
        double b = x + 1.0 - a, c = 1.0 / 1e-300, d = 1.0 / b, h = d;
        for (int i = 1; i <= 500; i++) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (fabs(d) < 1e-300) d = 1e-300;
            c = b + an / c;
            if (fabs(c) < 1e-300) c = 1e-300;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (fabs(del - 1.0) < 1e-14) break;
        }
        return exp(-x + a * log(x) - lgamma(a)) * h;
    }
}

/* Exact time-dependent Schrödinger propagator for a coherent double-slit wave packet */
typedef struct {
    double d;
    double sigma0;
    double T;
    double norm;
} wave_params_t;

static double complex eval_psi_inc(double x, const wave_params_t *wp) {
    double complex zt = 1.0 + I * (wp->T / (2.0 * wp->sigma0 * wp->sigma0));
    double pref = pow(2.0 * M_PI * wp->sigma0 * wp->sigma0, -0.25) / sqrt(wp->norm);
    double complex inv_z = 1.0 / zt;
    double complex csqrt_z = csqrt(zt);
    double x1 = x + 0.5 * wp->d;
    double x2 = x - 0.5 * wp->d;
    double complex g1 = cexp(-(x1 * x1) * inv_z / (4.0 * wp->sigma0 * wp->sigma0));
    double complex g2 = cexp(-(x2 * x2) * inv_z / (4.0 * wp->sigma0 * wp->sigma0));
    return (pref / csqrt_z) * (g1 + g2);
}

int main(void) {
    const int M = 60;                 /* Row of 60 solid-state detectors */
    const int N_SURF = 64;            /* Near-surface bound atomic sites per detector */
    const int K_CYC = 16;             /* Optical/de Broglie cycles during packet passage */
    const int N_TRIALS = 50000;       /* Number of single-photon arrivals */
    const double X_MIN = -15.0;
    const double X_MAX = +15.0;
    const double dx = (X_MAX - X_MIN) / (double)M;

    /* Wave packet parameters producing 7 interference fringes across [-15, +15] */
    wave_params_t wp;
    wp.d = 3.20;
    wp.sigma0 = 0.30;
    wp.T = 2.60;
    wp.norm = 2.0 * (1.0 + exp(-(wp.d * wp.d) / (8.0 * wp.sigma0 * wp.sigma0)));

    /* 8-point Gauss-Legendre nodes/weights on [-1, 1] */
    static const double gl_x[8] = {
        -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363
    };
    static const double gl_w[8] = {
         0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
         0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763
    };

    double x_center[M];
    double complex psi_det[M];
    double u_born[M];
    double cdf_born[M];
    double total_u = 0.0;

    for (int m = 0; m < M; m++) {
        double a = X_MIN + m * dx;
        double b = a + dx;
        double mid = 0.5 * (a + b);
        double half = 0.5 * (b - a);
        x_center[m] = mid;
        psi_det[m] = eval_psi_inc(mid, &wp);

        double integ = 0.0;
        for (int q = 0; q < 8; q++) {
            double xq = mid + half * gl_x[q];
            double complex pq = eval_psi_inc(xq, &wp);
            double dens = creal(pq) * creal(pq) + cimag(pq) * cimag(pq);
            integ += gl_w[q] * dens;
        }
        u_born[m] = half * integ;
        total_u += u_born[m];
    }

    double running = 0.0;
    for (int m = 0; m < M; m++) {
        u_born[m] /= total_u; /* Normalized incident energy share across detector row */
        running += u_born[m];
        cdf_born[m] = running;
    }
    cdf_born[M - 1] = 1.0;

    int *counts_copenhagen = calloc(M, sizeof(int));
    int *counts_tunneling  = calloc(M, sizeof(int));
    double *sum_dep_energy = calloc(M, sizeof(double));
    double *sum_noise_bias = calloc(M, sizeof(double));

    rng_t rng;
    rng_seed(&rng, 20261121ULL);

    /*
     * Surface bound-electron and vacuum fluctuation amplitudes (per cycle):
     * Weak zero-mean local background fluctuations around resonant frequency.
     */
    const double A_el = 0.008;
    const double sigma_vac = 0.008;

    /* Precompute 65536-point uniform unit-circle phasor table & carrier phases */
    static double complex phasor_lut[65536];
    for (int k = 0; k < 65536; k++) {
        double ang = 2.0 * M_PI * ((double)k + 0.5) / 65536.0;
        phasor_lut[k] = cos(ang) + I * sin(ang);
    }
    double complex carrier_lut[16];
    for (int c = 0; c < K_CYC; c++) {
        double ang = 2.0 * M_PI * (double)c / (double)K_CYC;
        carrier_lut[c] = cos(ang) - I * sin(ang);
    }

    double amp_site[60];
    double complex psi_site_arr[60];
    double u_site_arr[60];
    for (int m = 0; m < M; m++) {
        amp_site[m] = sqrt(u_born[m] / (double)N_SURF);
        double phi_inc = carg(psi_det[m]);
        psi_site_arr[m] = amp_site[m] * (cos(phi_inc) + I * sin(phi_inc));
        u_site_arr[m] = amp_site[m] * amp_site[m];
    }
    const double scale_el  = A_el / sqrt((double)N_SURF);
    const double scale_vac = sigma_vac / sqrt((double)N_SURF);

    /* Cascade phase-coherence tracking across K = 0, 1, 2, 5, 10, 50, 200 scatterings */
    const int K_STEPS[7] = {0, 1, 2, 5, 10, 50, 200};
    double complex coh_sum[7] = {0};
    int double_clicks = 0;

    for (int trial = 0; trial < N_TRIALS; trial++) {
        /* ------------------------------------------------------------------
         * PROTOCOL A: Copenhagen / Born-Rule + Direct RNG Selection
         * ------------------------------------------------------------------ */
        double u_rng = rng_uniform_open(&rng);
        int m_cop = M - 1;
        for (int m = 0; m < M; m++) {
            if (u_rng <= cdf_born[m]) {
                m_cop = m;
                break;
            }
        }
        counts_copenhagen[m_cop]++;

        /* ------------------------------------------------------------------
         * PROTOCOL B: Real-Wave Surface Interference + First-to-Tunnel Race
         *             + Nonlocal Photon Depletion + Solid-State Cascade
         * ------------------------------------------------------------------ */
        double best_time = 1e300;
        int win_m = -1;

        for (int m = 0; m < M; m++) {
            double complex psi_inc_site = psi_site_arr[m];
            double u_inc_site = u_site_arr[m];

            for (int a = 0; a < N_SURF; a++) {
                /* Pre-arrival sub-barrier action deficit (independent of psi_inc) */
                double delta_S = rng_uniform_open(&rng);

                /* Cycle-integrated deposited energy from 3-wave complex interference */
                double dep_accum = 0.0;
                double noise_accum = 0.0;
                for (int c = 0; c < K_CYC; c++) {
                    double complex w_inc = psi_inc_site * carrier_lut[c];

                    /* Draw 4 independent 16-bit random phase indices in one u64 call */
                    uint64_t r64 = rng_next_u64(&rng);
                    double complex w_el  = scale_el * phasor_lut[r64 & 0xFFFFU];
                    /* Two-mode zero-point vacuum superposition with random phases */
                    double complex w_vac = (0.7071067811865475 * scale_vac) *
                                           (phasor_lut[(r64 >> 16) & 0xFFFFU] +
                                            phasor_lut[(r64 >> 32) & 0xFFFFU]);

                    double complex w_bg  = w_el + w_vac;
                    double complex w_tot = w_inc + w_bg;

                    double e_tot = creal(w_tot) * creal(w_tot) + cimag(w_tot) * cimag(w_tot);
                    double e_bg  = creal(w_bg)  * creal(w_bg)  + cimag(w_bg)  * cimag(w_bg);

                    /* Net deposited energy above unilluminated bound+vacuum equilibrium */
                    double delta_u = e_tot - e_bg;
                    dep_accum += delta_u;
                    noise_accum += (delta_u - u_inc_site);
                }

                double u_bar = dep_accum / (double)K_CYC;
                sum_dep_energy[m] += u_bar;
                sum_noise_bias[m] += (noise_accum / (double)K_CYC);

                /* Only positive net deposited energy drives upward barrier tunneling */
                if (u_bar > 0.0) {
                    double t_tunnel = delta_S / u_bar;
                    if (t_tunnel < best_time) {
                        best_time = t_tunnel;
                        win_m = m;
                    } else if (t_tunnel == best_time) {
                        double_clicks++;
                    }
                }
            }
        }

        /* Nonlocal whole-photon absorption at winning tunneling site win_m */
        if (win_m >= 0) {
            counts_tunneling[win_m]++;

            /* Post-tunneling cascade of interactions with randomly phased solid-state electrons */
            double phi_packet = carg(psi_det[win_m]);
            double scrambled_phase = phi_packet;
            int k_curr = 0;
            for (int s = 0; s < 7; s++) {
                while (k_curr < K_STEPS[s]) {
                    double th_bound = 2.0 * M_PI * rng_uniform_open(&rng);
                    scrambled_phase += th_bound;
                    k_curr++;
                }
                /* Phase-retrieval correlation against original incident phase */
                coh_sum[s] += cexp(I * (scrambled_phase - phi_packet));
            }
        }
    }

    /* Compute Pearson chi-square and L1 metrics for both protocols */
    double chi2_cop = 0.0, chi2_tun = 0.0;
    double l1_cop = 0.0, l1_tun = 0.0;
    int dof = 0;

    FILE *fp = fopen("exp_v3_tunneling_vs_copenhagen.csv", "w");
    if (fp) {
        fprintf(fp, "pixel,x_center,born_prob,copenhagen_prob,tunneling_prob,mean_dep_energy,mean_noise_bias\n");
    }

    for (int m = 0; m < M; m++) {
        double p_b = u_born[m];
        double p_c = (double)counts_copenhagen[m] / (double)N_TRIALS;
        double p_t = (double)counts_tunneling[m]  / (double)N_TRIALS;
        double expected = p_b * (double)N_TRIALS;

        l1_cop += fabs(p_c - p_b);
        l1_tun += fabs(p_t - p_b);

        if (expected >= 5.0) {
            double dc = (double)counts_copenhagen[m] - expected;
            double dt = (double)counts_tunneling[m]  - expected;
            chi2_cop += (dc * dc) / expected;
            chi2_tun += (dt * dt) / expected;
            dof++;
        }

        if (fp) {
            fprintf(fp, "%d,%.4f,%.6f,%.6f,%.6f,%.6f,%.3e\n",
                    m, x_center[m], p_b, p_c, p_t,
                    sum_dep_energy[m] / (double)N_TRIALS,
                    sum_noise_bias[m] / (double)N_TRIALS);
        }
    }
    if (fp) fclose(fp);
    dof -= 1;

    double pval_cop = gamma_q(0.5 * dof, 0.5 * chi2_cop);
    double pval_tun = gamma_q(0.5 * dof, 0.5 * chi2_tun);

    printf("=== Born Selection v3.0 Verification (N = %d single-photon events, M = %d detectors) ===\n",
           N_TRIALS, M);
    printf("Protocol A (Copenhagen Born+RNG): chi2 = %.4f (dof = %d), p = %.4f, L1 = %.6f\n",
           chi2_cop, dof, pval_cop, l1_cop);
    printf("Protocol B (Tunneling+Cascade)  : chi2 = %.4f (dof = %d), p = %.4f, L1 = %.6f, double_clicks = %d\n\n",
           chi2_tun, dof, pval_tun, l1_tun, double_clicks);

    printf("Representative Fringe & Trough Detectors:\n");
    printf("  Det   x_center   Born_Exact   Copenhagen_RNG   Tunnel_Cascade   Mean_Noise_Bias\n");
    const int rep_bins[11] = {4, 9, 14, 19, 24, 29, 30, 35, 40, 45, 55};
    for (int i = 0; i < 11; i++) {
        int m = rep_bins[i];
        printf("  %3d   %+7.2f    %9.6f      %9.6f        %9.6f        %+.3e\n",
               m, x_center[m], u_born[m],
               (double)counts_copenhagen[m] / (double)N_TRIALS,
               (double)counts_tunneling[m] / (double)N_TRIALS,
               sum_noise_bias[m] / (double)N_TRIALS);
    }

    printf("\nCascade Phase-Information Irreversibility (Coherence Visibility |<exp(i*Delta_Phi)>|):\n");
    for (int s = 0; s < 7; s++) {
        double vis = cabs(coh_sum[s]) / (double)N_TRIALS;
        printf("  After K = %3d random-phase electron collisions: visibility = %.6f\n",
               K_STEPS[s], vis);
    }

    free(counts_copenhagen);
    free(counts_tunneling);
    free(sum_dep_energy);
    free(sum_noise_bias);
    return 0;
}
