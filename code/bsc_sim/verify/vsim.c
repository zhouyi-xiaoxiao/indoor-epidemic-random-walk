/* vsim.c -- second exact (direct-method) simulator written for the
 * re-check.  Deliberately different in structure from the
 * simulator gillespie.c:
 *   - the infection event is chosen as an (infective, kernel-cell) pair, so
 *     the infector is known directly;
 *   - removed people are not moved in mode 0 (they are irrelevant there);
 *   - a different RNG (splitmix64 stream);
 *   - all initial positions are drawn here (uniform over walkable cells)
 *     unless ipos >= 0.
 * mode 0: pair hazard kw(x,y) = beta q(x) P(y|x) / rho_bar
 * mode 1: per-cell rule, same cell, per-susceptible rate beta q(y) nI(y)/ntot(y)
 */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static uint64_t sm_state;
static inline uint64_t sm_next(void) {
    uint64_t z = (sm_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline double U(void) { return (sm_next() >> 11) * (1.0 / 9007199254740992.0); }
static inline int RI(int n) { int k = (int)(U() * n); return k < n ? k : n - 1; }

int vrun(int M, const int *nbr, const double *w,
         const int *kptr, const int *kidx, const double *kw, const double *bq,
         int mode, double gamma, int N, int n_rep, const int *ipos,
         const int *n_index_p, int gmax, double t_max, uint64_t seed,
         int nseg, const double *seg_end, const double *seg_hop, const double *seg_m,
         int *final, int *off, int *peak, double *peak_t, double *dur, double *expo,
         int *gens, int *cellgen, double *curve, double dt_curve, int n_curve,
         int move_all, int *n_unfinished)
{
    int n_index = n_index_p[0];
    double *wtot = (double *)calloc(M, sizeof(double));
    for (int x = 0; x < M; x++)
        for (int j = 0; j < 4; j++)
            if (nbr[4 * x + j] >= 0) wtot[x] += w[4 * x + j];
    int *pos = (int *)malloc(sizeof(int) * N);
    int *st = (int *)malloc(sizeof(int) * N);
    int *gen = (int *)malloc(sizeof(int) * N);
    int *nS = (int *)malloc(sizeof(int) * M);
    int *nI = (int *)malloc(sizeof(int) * M);
    int *nT = (int *)malloc(sizeof(int) * M);
    double period = seg_end[nseg - 1];
    int unfinished = 0;
    sm_state = seed;

    for (int rep = 0; rep < n_rep; rep++) {
        memset(nS, 0, sizeof(int) * M);
        memset(nI, 0, sizeof(int) * M);
        memset(nT, 0, sizeof(int) * M);
        for (int i = 0; i < N; i++) {
            int p = ipos[(size_t)rep * N + i];
            pos[i] = p >= 0 ? p : RI(M);
            st[i] = 0; gen[i] = -1;
            nT[pos[i]]++;
        }
        for (int i = 0; i < n_index; i++) { st[i] = 1; gen[i] = 0; }
        int cS = 0, cI = 0;
        for (int i = 0; i < N; i++) {
            if (st[i] == 0) { nS[pos[i]]++; cS++; } else { nI[pos[i]]++; cI++; }
        }
        int nfinal = n_index, noff = 0, pk = cI;
        double pkt = 0.0, t = 0.0, ex = 0.0;
        int seg = 0; double cyc = 0.0;
        for (int g = 0; g < 8; g++) gens[rep * 8 + g] = 0;
        gens[rep * 8] = n_index;
        int ck = 0;
        double text = NAN;

        while (1) {
            if (cI == 0) { text = t; break; }
            double hs = seg_hop[seg], m = seg_m[seg];
            /* --- rates, recomputed from scratch every event (brute force) */
            double H = 0.0;
            for (int i = 0; i < N; i++)
                if (move_all || st[i] != 2) H += wtot[pos[i]];
            H *= hs;
            double A = 0.0;
            if (mode == 0) {
                for (int i = 0; i < N; i++) if (st[i] == 1) {
                    int x = pos[i];
                    for (int k = kptr[x]; k < kptr[x + 1]; k++) A += kw[k] * nS[kidx[k]];
                }
            } else {
                for (int y = 0; y < M; y++)
                    if (nS[y] > 0 && nI[y] > 0) A += bq[y] * nS[y] * (double)nI[y] / nT[y];
            }
            A *= m;
            double tot = H + gamma * cI + A;
            double dt = -log(1.0 - U()) / tot;
            double tb = (nseg > 1) ? cyc + seg_end[seg] : INFINITY;
            double tn = t + dt;
            int brk = 0;
            if (tb < tn) { tn = tb; brk = 1; }
            if (t_max < tn) { tn = t_max; brk = 2; }
            while (ck < n_curve && ck * dt_curve < tn) { curve[(size_t)rep * n_curve + ck] = cI; ck++; }
            /* exposure of person 0 (potential rate if everybody in reach were susceptible) */
            if (st[0] == 1 && mode == 0) {
                int x = pos[0]; double r = 0.0;
                for (int k = kptr[x]; k < kptr[x + 1]; k++)
                    r += kw[k] * (nT[kidx[k]] - (kidx[k] == x ? 1 : 0));
                ex += m * r * (tn - t);
            }
            t = tn;
            if (brk == 2) { unfinished++; break; }
            if (brk == 1) { seg++; if (seg == nseg) { seg = 0; cyc += period; } continue; }
            double u = U() * tot;
            if (u < H) {                       /* hop */
                double acc = 0.0; int i = -1; double r = u / hs;
                for (int c = 0; c < N; c++) {
                    if (!move_all && st[c] == 2) continue;
                    acc += wtot[pos[c]]; i = c;
                    if (acc > r) break;
                }
                int x = pos[i];
                double r2 = U() * wtot[x], a2 = 0.0; int y = -1;
                for (int j = 0; j < 4; j++) {
                    if (nbr[4 * x + j] < 0 || w[4 * x + j] <= 0.0) continue;
                    y = nbr[4 * x + j]; a2 += w[4 * x + j];
                    if (a2 > r2) break;
                }
                if (y < 0) continue;
                pos[i] = y; nT[x]--; nT[y]++;
                if (st[i] == 0) { nS[x]--; nS[y]++; }
                else if (st[i] == 1) { nI[x]--; nI[y]++; }
            } else if (u < H + gamma * cI) {   /* recovery of a uniformly chosen infective */
                int k = RI(cI), i = -1;
                for (int c = 0; c < N; c++) if (st[c] == 1) { if (k == 0) { i = c; break; } k--; }
                st[i] = 2; nI[pos[i]]--; cI--;
            } else {                           /* infection */
                double r = (u - H - gamma * cI) / m, acc = 0.0;
                int src = -1, y = -1;
                if (mode == 0) {
                    for (int i = 0; i < N && acc <= r; i++) if (st[i] == 1) {
                        int x = pos[i];
                        for (int k = kptr[x]; k < kptr[x + 1]; k++) {
                            double a = kw[k] * nS[kidx[k]];
                            if (a <= 0.0) continue;
                            acc += a; src = i; y = kidx[k];
                            if (acc > r) break;
                        }
                    }
                } else {
                    for (int c = 0; c < M; c++) {
                        if (nS[c] > 0 && nI[c] > 0) {
                            acc += bq[c] * nS[c] * (double)nI[c] / nT[c]; y = c;
                            if (acc > r) break;
                        }
                    }
                    int k = RI(nI[y]);
                    for (int c = 0; c < N; c++) if (st[c] == 1 && pos[c] == y) { if (k == 0) { src = c; break; } k--; }
                }
                int k = RI(nS[y]), j = -1;
                for (int c = 0; c < N; c++) if (st[c] == 0 && pos[c] == y) { if (k == 0) { j = c; break; } k--; }
                int g = gen[src] + 1;
                gen[j] = g; nfinal++; cS--; nS[y]--;
                if (src == 0) noff++;
                if (g < 8) gens[rep * 8 + g]++;
                if (g >= 1 && g <= 4) cellgen[(g - 1) * M + y]++;
                if (gmax >= 0 && g > gmax) { st[j] = 2; }
                else { st[j] = 1; nI[y]++; cI++; if (cI > pk) { pk = cI; pkt = t; } }
            }
        }
        while (ck < n_curve) { curve[(size_t)rep * n_curve + ck] = cI; ck++; }
        final[rep] = nfinal; off[rep] = noff; peak[rep] = pk; peak_t[rep] = pkt;
        dur[rep] = text; expo[rep] = ex;
    }
    n_unfinished[0] = unfinished;
    free(wtot); free(pos); free(st); free(gen); free(nS); free(nI); free(nT);
    return 0;
}
