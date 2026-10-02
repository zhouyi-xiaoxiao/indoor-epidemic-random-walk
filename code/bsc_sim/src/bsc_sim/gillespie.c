/*
 * gillespie.c -- exact event-driven (Gillespie / direct method) simulation of
 * the lattice random-walk SIR model with explicit individuals.
 *
 * Events
 *   hop        person i at cell x jumps to neighbour y      rate W(x -> y)
 *   recovery   infective -> removed                         rate gamma
 *   infection  susceptible at cell y becomes infective      rate m(t) * E(y)
 *
 * Transmission rule (mode 0, pair hazard): every infective at cell x
 * adds c(x, y) = beta q(x) P(y|x) / rho_bar to the exposure field E(y) of
 * each cell y in its contact neighbourhood; every susceptible at y is
 * infected at rate m(t) E(y).  P is a row-stochastic contact kernel and
 * rho_bar the mean number of other people per cell, so that an infective in
 * a fully susceptible room infects at the expected rate beta q(x).
 *
 * Transmission rule (mode 1, per-cell rule): same cell only, rate per
 * susceptible beta q(y) n_I(y) / n_total(y).  Used for the remark that this
 * rule is sub-critical at the occupancies of the scenes.
 *
 * Time dependence: a periodic, piecewise-constant schedule selects a hop-rate
 * profile and a transmission multiplier m for each segment.  Because all
 * rates are constant inside a segment and the process is Markov, restarting
 * the clock at a segment boundary is exact.
 *
 * Build:  clang -O2 -shared -fPIC -o libgillespie.dylib gillespie.c
 */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------- RNG */
typedef struct { uint64_t s[4]; } rng_t;

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}
static inline uint64_t rng_next(rng_t *r) {           /* xoshiro256** */
    const uint64_t result = rotl(r->s[1] * 5, 7) * 9;
    const uint64_t t = r->s[1] << 17;
    r->s[2] ^= r->s[0];
    r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2];
    r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3] = rotl(r->s[3], 45);
    return result;
}
static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static void rng_seed(rng_t *r, uint64_t seed) {
    for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&seed);
}
static inline double rng_u(rng_t *r) {                 /* uniform on [0,1) */
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}
static inline int rng_int(rng_t *r, int n) {           /* uniform on 0..n-1 */
    int k = (int)(rng_u(r) * n);
    return k < n ? k : n - 1;
}

/* ---------------------------------------------------------------- state */
typedef struct {
    int M, N, mode;
    const int *ker_ptr, *ker_idx;
    const double *ker_w, *bq;
    int *pos, *nS, *nI, *ntot, *inf_list, n_inf;
    signed char *state;       /* 0 S, 1 I (infectious), 2 removed */
    double *E, *a, A;
} sim_t;

static inline void upd_site(sim_t *s, int y) {
    double v;
    if (s->mode == 0) {
        v = s->nS[y] * s->E[y];
        if (v < 0.0) v = 0.0;
    } else {
        v = (s->nS[y] > 0 && s->nI[y] > 0)
                ? s->bq[y] * s->nS[y] * s->nI[y] / (double)s->ntot[y] : 0.0;
    }
    s->A += v - s->a[y];
    s->a[y] = v;
}

static void recompute_inf(sim_t *s) {
    if (s->mode == 0) {
        memset(s->E, 0, sizeof(double) * s->M);
        for (int k = 0; k < s->n_inf; k++) {
            int x = s->pos[s->inf_list[k]];
            for (int j = s->ker_ptr[x]; j < s->ker_ptr[x + 1]; j++)
                s->E[s->ker_idx[j]] += s->ker_w[j];
        }
    }
    s->A = 0.0;
    for (int y = 0; y < s->M; y++) {
        double v;
        if (s->mode == 0)
            v = s->nS[y] * s->E[y];
        else
            v = (s->nS[y] > 0 && s->nI[y] > 0)
                    ? s->bq[y] * s->nS[y] * s->nI[y] / (double)s->ntot[y] : 0.0;
        s->a[y] = v;
        s->A += v;
    }
}

/* potential (mean-field) infection rate of one infective: the rate at which
 * it would infect if everybody within its contact neighbourhood were
 * susceptible.  Its time integral has expectation (1^T K)_x exactly. */
static double potential_rate(const sim_t *s, int agent) {
    int x = s->pos[agent];
    if (s->mode == 0) {
        double r = 0.0;
        for (int j = s->ker_ptr[x]; j < s->ker_ptr[x + 1]; j++) {
            int y = s->ker_idx[j];
            r += s->ker_w[j] * (s->ntot[y] - (y == x ? 1 : 0));
        }
        return r;
    }
    return s->bq[x] * (s->ntot[x] - 1) / (double)s->ntot[x];
}

/* ---------------------------------------------------------------- driver */
int run_batch(
    int M, const int *nbr_ptr, const int *nbr_idx,
    int n_prof, const double *W,
    const int *ker_ptr, const int *ker_idx, const double *ker_w,
    const double *bq, int mode, double gamma,
    int n_seg, const double *seg_end, const int *seg_prof, const double *seg_m,
    int N, int n_rep, const int *init_pos, int n_index, const int *index_agent,
    int gmax, double t_max, double dt_out, int n_out, uint64_t seed,
    int n_snap, const int *snap_idx,
    int32_t *out_S, int32_t *out_I,
    float *t_inf, float *t_rec, int32_t *site_inf, int32_t *gen,
    int32_t *infector, double *index_exposure, double *t_ext,
    uint16_t *snap, int64_t *n_events)
{
    const int nnz = nbr_ptr[M];
    double *Wtot = (double *)calloc((size_t)n_prof * M, sizeof(double));
    double *Wmax = (double *)calloc((size_t)n_prof, sizeof(double));
    for (int p = 0; p < n_prof; p++)
        for (int x = 0; x < M; x++) {
            double t = 0.0;
            for (int j = nbr_ptr[x]; j < nbr_ptr[x + 1]; j++)
                t += W[(size_t)p * nnz + j];
            Wtot[(size_t)p * M + x] = t;
            if (t > Wmax[p]) Wmax[p] = t;
        }

    sim_t s;
    s.M = M; s.N = N; s.mode = mode;
    s.ker_ptr = ker_ptr; s.ker_idx = ker_idx; s.ker_w = ker_w; s.bq = bq;
    s.pos = (int *)malloc(sizeof(int) * N);
    s.state = (signed char *)malloc(N);
    s.inf_list = (int *)malloc(sizeof(int) * N);
    s.nS = (int *)malloc(sizeof(int) * M);
    s.nI = (int *)malloc(sizeof(int) * M);
    s.ntot = (int *)malloc(sizeof(int) * M);
    s.E = (double *)malloc(sizeof(double) * M);
    s.a = (double *)malloc(sizeof(double) * M);
    const double period = seg_end[n_seg - 1];
    const double INF = 1e300;

    for (int rep = 0; rep < n_rep; rep++) {
        rng_t rng;
        rng_seed(&rng, seed + 0x632BE59BD9B4E019ULL * (uint64_t)(rep + 1));
        float *ti = t_inf + (size_t)rep * N, *tr = t_rec + (size_t)rep * N;
        int32_t *si = site_inf + (size_t)rep * N, *ge = gen + (size_t)rep * N;
        int32_t *io = infector + (size_t)rep * N;
        int32_t *oS = out_S + (size_t)rep * n_out, *oI = out_I + (size_t)rep * n_out;

        memset(s.nS, 0, sizeof(int) * M);
        memset(s.nI, 0, sizeof(int) * M);
        memset(s.ntot, 0, sizeof(int) * M);
        for (int i = 0; i < N; i++) {
            s.pos[i] = init_pos[(size_t)rep * N + i];
            s.state[i] = 0;
            s.nS[s.pos[i]]++;
            s.ntot[s.pos[i]]++;
            ti[i] = -1.0f; tr[i] = -1.0f; si[i] = -1; ge[i] = -1; io[i] = -1;
        }
        s.n_inf = 0;
        int nS_tot = N, n_removed = 0;
        for (int k = 0; k < n_index; k++) {
            int a = index_agent[(size_t)rep * n_index + k];
            if (s.state[a] != 0) continue;
            s.state[a] = 1;
            s.nS[s.pos[a]]--; s.nI[s.pos[a]]++;
            s.inf_list[s.n_inf++] = a;
            nS_tot--;
            ti[a] = 0.0f; ge[a] = 0; si[a] = s.pos[a];
        }
        const int idx0 = index_agent[(size_t)rep * n_index];
        double expo = 0.0;

        int seg = 0, prof = seg_prof[0];
        double m = seg_m[0], cyc = 0.0;
        double next_break = (n_seg > 1) ? seg_end[0] : INF;
        const double *Wp = W + (size_t)prof * nnz;
        const double *Wt = Wtot + (size_t)prof * M;
        double H = 0.0;
        for (int i = 0; i < N; i++) H += Wt[s.pos[i]];
        recompute_inf(&s);

        double t = 0.0;
        int out_k = 0, snap_k = 0;
        int64_t nev = 0, hop_since = 0;
        double text = NAN;

        while (1) {
            if (s.n_inf == 0) { text = t; break; }
            double total = H + gamma * s.n_inf + m * s.A;
            double t_new = t - log(1.0 - rng_u(&rng)) / total;
            double t_stop = t_new;
            int kind = 0;                     /* 0 event, 1 break, 2 end */
            if (next_break < t_stop) { t_stop = next_break; kind = 1; }
            if (t_max <= t_stop) { t_stop = t_max; kind = 2; }
            /* record the state that holds on [t, t_stop) */
            while (out_k < n_out && out_k * dt_out < t_stop) {
                oS[out_k] = nS_tot; oI[out_k] = s.n_inf;
                if (snap_k < n_snap && snap_idx[snap_k] == out_k) {
                    uint16_t *sn = snap + ((size_t)rep * n_snap + snap_k) * M;
                    for (int y = 0; y < M; y++) sn[y] = (uint16_t)s.nI[y];
                    snap_k++;
                }
                out_k++;
            }
            if (s.state[idx0] == 1)
                expo += m * potential_rate(&s, idx0) * (t_stop - t);
            t = t_stop;
            if (kind == 2) break;
            if (kind == 1) {
                seg++;
                if (seg == n_seg) { seg = 0; cyc += period; }
                next_break = cyc + seg_end[seg];
                m = seg_m[seg];
                if (seg_prof[seg] != prof) {
                    prof = seg_prof[seg];
                    Wp = W + (size_t)prof * nnz;
                    Wt = Wtot + (size_t)prof * M;
                    H = 0.0;
                    for (int i = 0; i < N; i++) H += Wt[s.pos[i]];
                }
                continue;
            }
            nev++;
            double u = rng_u(&rng) * total;
            if (u < H) {
                /* ---- hop: pick a person with probability ~ exit rate */
                int i = -1;
                for (int tries = 0; tries < 256; tries++) {
                    int c = rng_int(&rng, N);
                    if (rng_u(&rng) * Wmax[prof] < Wt[s.pos[c]]) { i = c; break; }
                }
                if (i < 0) {                   /* rare: exact linear scan */
                    double r = rng_u(&rng) * H, acc = 0.0;
                    i = N - 1;
                    for (int c = 0; c < N; c++) {
                        acc += Wt[s.pos[c]];
                        if (acc > r) { i = c; break; }
                    }
                }
                int x = s.pos[i];
                double r = rng_u(&rng) * Wt[x], acc = 0.0;
                int x2 = -1;
                for (int j = nbr_ptr[x]; j < nbr_ptr[x + 1]; j++) {
                    if (Wp[j] <= 0.0) continue;
                    x2 = nbr_idx[j];
                    acc += Wp[j];
                    if (acc > r) break;
                }
                if (x2 < 0) continue;
                s.pos[i] = x2;
                s.ntot[x]--; s.ntot[x2]++;
                H += Wt[x2] - Wt[x];
                if (s.state[i] == 0) {
                    s.nS[x]--; s.nS[x2]++;
                    upd_site(&s, x); upd_site(&s, x2);
                } else if (s.state[i] == 1) {
                    s.nI[x]--; s.nI[x2]++;
                    if (mode == 0) {
                        for (int j = ker_ptr[x]; j < ker_ptr[x + 1]; j++) {
                            s.E[ker_idx[j]] -= ker_w[j];
                            upd_site(&s, ker_idx[j]);
                        }
                        for (int j = ker_ptr[x2]; j < ker_ptr[x2 + 1]; j++) {
                            s.E[ker_idx[j]] += ker_w[j];
                            upd_site(&s, ker_idx[j]);
                        }
                    } else {
                        upd_site(&s, x); upd_site(&s, x2);
                    }
                } else if (mode == 1) {
                    upd_site(&s, x); upd_site(&s, x2);
                }
                if (++hop_since >= 20000) {    /* bound floating-point drift */
                    hop_since = 0;
                    recompute_inf(&s);
                    H = 0.0;
                    for (int c = 0; c < N; c++) H += Wt[s.pos[c]];
                }
            } else if (u < H + gamma * s.n_inf) {
                /* ---- recovery */
                int k = rng_int(&rng, s.n_inf);
                int j = s.inf_list[k];
                s.inf_list[k] = s.inf_list[--s.n_inf];
                s.state[j] = 2;
                s.nI[s.pos[j]]--;
                tr[j] = (float)t;
                n_removed++;
                recompute_inf(&s);
            } else {
                /* ---- infection: cell ~ a[y], then a susceptible there */
                double r = rng_u(&rng) * s.A, acc = 0.0;
                int y = -1;
                for (int c = 0; c < M; c++) {
                    if (s.a[c] <= 0.0) continue;
                    y = c;
                    acc += s.a[c];
                    if (acc > r) break;
                }
                if (y < 0 || s.nS[y] <= 0) { recompute_inf(&s); continue; }
                int k = rng_int(&rng, s.nS[y]), j = -1;
                for (int c = 0; c < N; c++)
                    if (s.pos[c] == y && s.state[c] == 0) {
                        if (k == 0) { j = c; break; }
                        k--;
                    }
                /* infector ~ its contribution to E(y) */
                int src = -1;
                if (mode == 0) {
                    double rr = rng_u(&rng) * s.E[y], ac = 0.0;
                    for (int c = 0; c < s.n_inf; c++) {
                        int ia = s.inf_list[c], x = s.pos[ia];
                        for (int q = ker_ptr[x]; q < ker_ptr[x + 1]; q++)
                            if (ker_idx[q] == y) {
                                ac += ker_w[q];
                                src = ia;
                                break;
                            }
                        if (ac > rr) break;
                    }
                } else {
                    int kk = rng_int(&rng, s.nI[y]);
                    for (int c = 0; c < s.n_inf; c++) {
                        int ia = s.inf_list[c];
                        if (s.pos[ia] == y) {
                            if (kk == 0) { src = ia; break; }
                            kk--;
                        }
                    }
                }
                int g = (src >= 0 ? ge[src] : 0) + 1;
                s.nS[y]--; nS_tot--;
                ti[j] = (float)t; si[j] = y; ge[j] = g; io[j] = src;
                if (gmax >= 0 && g > gmax) {   /* counted, but not infectious */
                    s.state[j] = 2;
                    tr[j] = (float)t;
                    n_removed++;
                } else {
                    s.state[j] = 1;
                    s.nI[y]++;
                    s.inf_list[s.n_inf++] = j;
                }
                recompute_inf(&s);
            }
        }
        while (out_k < n_out) {
            oS[out_k] = nS_tot; oI[out_k] = s.n_inf;
            if (snap_k < n_snap && snap_idx[snap_k] == out_k) {
                uint16_t *sn = snap + ((size_t)rep * n_snap + snap_k) * M;
                for (int y = 0; y < M; y++) sn[y] = (uint16_t)s.nI[y];
                snap_k++;
            }
            out_k++;
        }
        index_exposure[rep] = expo;
        t_ext[rep] = text;
        n_events[rep] = nev;
    }
    free(Wtot); free(Wmax);
    free(s.pos); free(s.state); free(s.inf_list);
    free(s.nS); free(s.nI); free(s.ntot); free(s.E); free(s.a);
    return 0;
}
