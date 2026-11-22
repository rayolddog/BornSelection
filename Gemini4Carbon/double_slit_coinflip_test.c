/*
 * double_slit_coinflip_test.c
 *
 * Tests whether fair 50/50 coin flips alone (driving the multi-pixel martingale
 * selection dynamics at the detector screen) can reproduce BOTH:
 *   1. Multiple orders of interference fringes (7 fringes: m = -3,-2,-1,0,+1,+2,+3),
 *   2. The continuous decaying diffraction envelope of peak heights (100%, 70%, 25%, 4%),
 *   3. The smooth continuous shape of every individual peak and dark trough.
 *
 * Compile & run:
 *   clang -O3 -Wall -Wextra -std=c11 double_slit_coinflip_test.c -lm -o double_slit_coinflip_test
 *   ./double_slit_coinflip_test
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

typedef struct {
    double d;
    double sigma0;
    double T;
    double norm_fac;
} dslit_t;

static double dslit_born_density(double x, const dslit_t *ctx) {
    const double s2 = ctx->sigma0 * ctx->sigma0;
    const double complex z = 1.0 + I * (ctx->T / (2.0 * s2));
    const double x1 = -0.5 * ctx->d;
    const double x2 = +0.5 * ctx->d;
    const double pref = pow(2.0 * M_PI * s2, -0.25) / sqrt(ctx->norm_fac);
    const double complex z_inv_sqrt = 1.0 / csqrt(z);

    double complex psi1 = pref * z_inv_sqrt * cexp(-((x - x1) * (x - x1)) / (4.0 * s2 * z));
    double complex psi2 = pref * z_inv_sqrt * cexp(-((x - x2) * (x - x2)) / (4.0 * s2 * z));
    double complex psi = psi1 + psi2;
    return creal(psi) * creal(psi) + cimag(psi) * cimag(psi);
}

static double quad_gl8_dslit(const dslit_t *ctx, double a, double b) {
    static const double xi[4] = {
        0.183434642495649805, 0.525532409916328986,
        0.796666477413626740, 0.960289856497536232
    };
    static const double wi[4] = {
        0.362683783378361983, 0.313706645877887287,
        0.222381034453374471, 0.101228536290376259
    };
    double mid = 0.5 * (a + b), half = 0.5 * (b - a), sum = 0.0;
    for (int i = 0; i < 4; i++) {
        double dx = half * xi[i];
        sum += wi[i] * (dslit_born_density(mid - dx, ctx) + dslit_born_density(mid + dx, ctx));
    }
    return half * sum;
}

#define N_PIXELS 80

int main(void) {
    rng_t rng;
    rng_seed(&rng, 20261117ULL);

    /* 7-fringe Young double-slit configuration:
       Fringe spacing ~ 2*pi*T/d = 5.39 -> peaks at 0, +-5.3, +-10.6, +-15.8
       Diffraction envelope width Sigma(T) ~ 6.43 -> peak heights 100%, 70%, 25%, 4.3% */
    const double d = 4.20;
    const double sigma0 = 0.28;
    const double T = 3.60;
    const double norm_fac = 2.0 * (1.0 + exp(-(d * d) / (8.0 * sigma0 * sigma0)));
    const dslit_t ctx = { d, sigma0, T, norm_fac };

    const double x_min = -20.0;
    const double x_max = +20.0;
    const double dx = (x_max - x_min) / (double)N_PIXELS; /* 0.50 per pixel */

    /* Wave propagation through both slits produces initial screen weights p0[m] */
    double p0[N_PIXELS];
    double sum_p0 = 0.0;
    for (int m = 0; m < N_PIXELS; m++) {
        double a = x_min + m * dx;
        p0[m] = quad_gl8_dslit(&ctx, a, a + dx);
        sum_p0 += p0[m];
    }
    /* Normalize across the 80 screen pixels (captures >99.8% of total wave) */
    for (int m = 0; m < N_PIXELS; m++) p0[m] /= sum_p0;

    const int n_photons = 200000;
    const double eps = 0.28;      /* Martingale coupling (< 0.5 guarantees positivity) */
    const int max_steps = 600;
    long counts[N_PIXELS] = {0};
    long unresolved = 0;

    for (int shot = 0; shot < n_photons; shot++) {
        double p[N_PIXELS];
        memcpy(p, p0, sizeof(p));

        int winner = -1;
        for (int step = 0; step < max_steps; step++) {
            /* Draw 80 independent FAIR 50/50 coin flips xi_m in {-1, +1} using 2 u64 words */
            uint64_t bits0 = rng_next_u64(&rng);
            uint64_t bits1 = rng_next_u64(&rng);
            double xi[N_PIXELS];
            double xi_bar = 0.0;

            for (int m = 0; m < 64; m++) {
                double s = ((bits0 >> m) & 1ULL) ? 1.0 : -1.0;
                xi[m] = s;
                xi_bar += p[m] * s;
            }
            for (int m = 64; m < N_PIXELS; m++) {
                double s = ((bits1 >> (m - 64)) & 1ULL) ? 1.0 : -1.0;
                xi[m] = s;
                xi_bar += p[m] * s;
            }

            double p_max = 0.0;
            int argmax = 0;
            for (int m = 0; m < N_PIXELS; m++) {
                p[m] *= (1.0 + eps * (xi[m] - xi_bar));
                if (p[m] > p_max) {
                    p_max = p[m];
                    argmax = m;
                }
            }
            if (p_max > 0.999999) {
                winner = argmax;
                break;
            }
        }
        if (winner >= 0) counts[winner]++;
        else unresolved++;
    }

    int resolved = n_photons - (int)unresolved;
    printf("7-Fringe Single-Photon Double-Slit via Fair 50/50 Coin Flips\n");
    printf("Photons: %d | Resolved single-pixel clicks: %d | Unresolved: %ld\n\n",
           n_photons, resolved, unresolved);

    /* Print fringe peak heights (m = -3, -2, -1, 0, +1, +2, +3) and dark troughs */
    const int key_bins[13] = {
        8,   /* x = -15.75 (m = -3 outer fringe peak) */
        13,  /* x = -13.25 (dark trough) */
        19,  /* x = -10.25 (m = -2 fringe peak) */
        24,  /* x =  -7.75 (dark trough) */
        29,  /* x =  -5.25 (m = -1 fringe peak) */
        34,  /* x =  -2.75 (dark trough) */
        40,  /* x =  +0.25 (m =  0 central peak) */
        45,  /* x =  +2.75 (dark trough) */
        50,  /* x =  +5.25 (m = +1 fringe peak) */
        55,  /* x =  +7.75 (dark trough) */
        60,  /* x = +10.25 (m = +2 fringe peak) */
        66,  /* x = +13.25 (dark trough) */
        71   /* x = +15.75 (m = +3 outer fringe peak) */
    };
    const char *labels[13] = {
        "m = -3 peak  ", "dark trough  ", "m = -2 peak  ", "dark trough  ",
        "m = -1 peak  ", "dark trough  ", "m =  0 center", "dark trough  ",
        "m = +1 peak  ", "dark trough  ", "m = +2 peak  ", "dark trough  ",
        "m = +3 peak  "
    };

    printf("Feature        Pixel x    Born Prob     Coin-Flip Prob   Clicks\n");
    printf("-------------  -------    ---------     --------------   ------\n");
    for (int k = 0; k < 13; k++) {
        int b = key_bins[k];
        double x_mid = x_min + (b + 0.5) * dx;
        double p_emp = (double)counts[b] / (double)resolved;
        printf("%s  %+6.2f     %9.6f       %9.6f      %6ld\n",
               labels[k], x_mid, p0[b], p_emp, counts[b]);
    }

    double chi2 = 0.0, l1 = 0.0;
    int active = 0;
    double merged_exp = 0.0;
    long merged_obs = 0;
    for (int m = 0; m < N_PIXELS; m++) {
        double exp_m = p0[m] * (double)resolved;
        double p_emp = (double)counts[m] / (double)resolved;
        l1 += fabs(p_emp - p0[m]);
        if (exp_m >= 20.0) {
            double d_cnt = (double)counts[m] - exp_m;
            chi2 += (d_cnt * d_cnt) / exp_m;
            active++;
        } else {
            merged_exp += exp_m;
            merged_obs += counts[m];
        }
    }
    if (merged_exp >= 5.0) {
        double d_cnt = (double)merged_obs - merged_exp;
        chi2 += (d_cnt * d_cnt) / merged_exp;
        active++;
    }
    int dof = active - 1;
    printf("\nAll 80 Pixels Chi-Square: chi^2 = %.4f (dof = %d), p-value = %.4f, L1 = %.6f\n",
           chi2, dof, gamma_q(0.5 * dof, 0.5 * chi2), l1);
    return 0;
}
