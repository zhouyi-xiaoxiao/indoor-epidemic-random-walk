"""verify_theorems.py -- numerical verification of the theorems on small lattices (Supplementary Table S40 of the article).

Run:  python verify_theorems.py            (writes ../results/verify_theorems.json and .log)
Every check reports the worst-case violation over all random instances; nothing is tuned.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).parent))
import ngm_lattice as ng

OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20261001)
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


# ------------------------------------------------------------------ random instances


def random_lattice(rng, nmin=3, nmax=8, max_obst=0.3):
    while True:
        ny, nx = rng.integers(nmin, nmax + 1, size=2)
        mask = rng.random((ny, nx)) > rng.uniform(0, max_obst)
        lat = ng.Lattice(mask)
        if lat.n < 4:
            continue
        ncomp, lab = connected_components(lat.adjacency(), directed=False)
        if ncomp > 1:  # keep the largest component
            big = np.argmax(np.bincount(lab))
            m2 = np.zeros_like(mask)
            cells = lat.cells[lab == big]
            m2[cells[:, 0], cells[:, 1]] = True
            lat = ng.Lattice(m2)
        if lat.n >= 4:
            return lat


def random_q(rng, n):
    kind = rng.integers(0, 4)
    if kind == 0:
        q = rng.lognormal(0, 0.6, n)
    elif kind == 1:  # zoned: a few levels
        levels = rng.uniform(0.3, 3.0, rng.integers(2, 5))
        q = levels[rng.integers(0, len(levels), n)]
    elif kind == 2:  # some zeros
        q = rng.uniform(0.2, 2.5, n) * (rng.random(n) > 0.3)
        if q.max() == 0:
            q[rng.integers(n)] = 1.0
    else:
        q = rng.uniform(0.5, 2.0, n)
    return q


def random_instance(rng, kind=None):
    kind = kind or rng.choice(["symmetric", "departure", "nonrev", "dense"])
    if kind == "dense":
        n = int(rng.integers(3, 12))
        W = rng.lognormal(0, 1, (n, n)) * (rng.random((n, n)) < 0.7)
        np.fill_diagonal(W, 0)
        # make irreducible by adding a cycle
        for i in range(n):
            W[(i + 1) % n, i] += rng.uniform(0.1, 1)
        M = W - np.diag(W.sum(axis=0))
        M = sp.csr_matrix(M)
        lat = None
    else:
        lat = random_lattice(rng)
        n = lat.n
        D = rng.lognormal(0, 0.8, n)
        if kind == "symmetric":
            M = ng.generator_symmetric(lat, D, kind=rng.choice(["harmonic", "arithmetic"]))
        elif kind == "departure":
            M = ng.generator_departure(lat, D)
        else:
            wf = rng.lognormal(0, 1, len(lat.edges))
            wb = rng.lognormal(0, 1, len(lat.edges))
            M = ng.generator_from_edge_rates(n, lat.edges, wf, wb)
    q = random_q(rng, n)
    gamma = float(rng.uniform(0.05, 2.0))
    beta = float(rng.uniform(0.02, 3.0))
    pi = ng.stationary(M)
    return dict(kind=kind, lat=lat, M=M, q=q, beta=beta, gamma=gamma, pi=pi, n=n)


RES = {}
t0 = time.time()

# ------------------------------------------------------------------ Theorem 1
log("== Theorem 1: NGM, sign equivalence, root characterisation ==")
N1 = 400
worst_root, sign_fail, worst_pi, worst_thr = 0.0, 0, 0.0, 0.0
worst_imag = 0.0
kinds_count = {}
T1_PTS = []
for it in range(N1):
    ins = random_instance(rng)
    kinds_count[ins["kind"]] = kinds_count.get(ins["kind"], 0) + 1
    M, q, b, g = ins["M"], ins["q"], ins["beta"], ins["gamma"]
    R0 = ng.R0_dense(M, q, b, g)
    lam, u, l, im = ng.growth_rate_dense(M, q, b, g)
    worst_imag = max(worst_imag, im)
    # (a) root characterisation: s(M + bQ/R0 - g) = 0
    worst_root = max(worst_root, abs(ng.root_characterisation(M, q, b, g, R0)))
    T1_PTS.append((R0, lam / g, ["symmetric", "departure", "nonrev", "dense"].index(str(ins["kind"]))))
    # (b) sign equivalence (skip numerically-critical cases)
    if abs(R0 - 1) > 1e-9 and np.sign(R0 - 1) != np.sign(lam):
        sign_fail += 1
    # (c) exactly at threshold: rescale beta so that R0 = 1  ->  lambda_1 = 0
    lam_thr, *_ = ng.growth_rate_dense(M, q, b / R0, g)
    worst_thr = max(worst_thr, abs(lam_thr))
    # stationary distribution sanity
    worst_pi = max(worst_pi, np.abs(M @ ins["pi"]).max())
np.save(OUT / "T1_points.npy", np.array(T1_PTS))
kinds_count = {str(k): v for k, v in kinds_count.items()}
log(f"instances {N1} by kind {kinds_count}")
log(f"max |s(M+bQ/R0-g)|           = {worst_root:.3e}")
log(f"sign(R0-1) != sign(lambda1)  : {sign_fail} failures")
log(f"max |lambda1| at beta/R0     = {worst_thr:.3e}")
log(f"max imag part of Perron root = {worst_imag:.3e}")
RES["T1"] = dict(n_instances=N1, kinds=kinds_count, max_abs_root_residual=worst_root,
                 sign_failures=sign_fail, max_abs_lambda_at_threshold=worst_thr)

# (d) K_rs = beta q_r * int_0^inf e^{-gamma t} P(r,t|s) dt : quadrature with expm
log("-- Theorem 1(d): K = beta Q * Laplace transform of the propagator at gamma")
ins = random_instance(rng, "symmetric")
M, q, b, g = ins["M"].toarray(), ins["q"], ins["beta"], ins["gamma"]
n = ins["n"]
K, Vinv = ng.ngm_dense(M, q, b, g)
# Gauss-Laguerre quadrature in tau = gamma t
xs, ws = np.polynomial.laguerre.laggauss(80)
Gint = sum(w * sla.expm(M * x / g) for x, w in zip(xs, ws)) / g
err_lap = np.abs(Gint - Vinv).max() / np.abs(Vinv).max()
log(f"n={n}: rel. error of Laguerre-quadrature Laplace transform vs (gamma-M)^-1 = {err_lap:.2e}")
# column sums of expm(Mt) are 1 (conservation); 1^T V^{-1} = 1^T/gamma
log(f"max |1^T V^-1 - 1/gamma| = {np.abs(Vinv.sum(axis=0) - 1 / g).max():.2e}")
RES["T1d_laplace_rel_err"] = float(err_lap)

# (e) Monte-Carlo check of the sojourn-time interpretation (single infective, no epidemic)
log("-- Theorem 1(e): Monte-Carlo sojourn times of one infective vs V^-1 column")
lat = ng.Lattice(np.ones((3, 4), bool))
Dm = rng.uniform(0.5, 2.0, lat.n)
Mmc = ng.generator_symmetric(lat, Dm).toarray()
gmc = 0.7
Vinv_mc = np.linalg.inv(gmc * np.eye(lat.n) - Mmc)
s0 = 5
nbatch, per = 40, 100000
nrep = nbatch * per
rates_out = -np.diag(Mmc)
Pjump = Mmc.copy()
np.fill_diagonal(Pjump, 0)
Pjump = Pjump / Pjump.sum(axis=0, keepdims=True)  # column-stochastic jump chain
cum = np.cumsum(Pjump, axis=0)
soj_b = np.zeros((nbatch, lat.n))
rr = np.random.default_rng(20261001)
for bidx in range(nbatch):  # vectorised Gillespie: all walkers of a batch advance one event at a time
    pos = np.full(per, s0)
    rem = rr.exponential(1 / gmc, per)  # remaining infectious lifetime
    alive = np.ones(per, bool)
    while alive.any():
        idx = np.nonzero(alive)[0]
        hold = rr.exponential(1 / rates_out[pos[idx]])
        dt = np.minimum(hold, rem[idx])
        np.add.at(soj_b[bidx], pos[idx], dt)
        rem[idx] -= dt
        dies = dt < hold
        alive[idx[dies]] = False
        mv = idx[~dies]
        pos[mv] = (rr.random(len(mv))[None, :] > cum[:, pos[mv]]).sum(axis=0)
    soj_b[bidx] /= per
soj = soj_b.mean(axis=0)
se = soj_b.std(axis=0, ddof=1) / np.sqrt(nbatch)
zsc = (soj - Vinv_mc[:, s0]) / se
zmax = np.abs(soj - Vinv_mc[:, s0]).max()
log(f"3x4 lattice, {nrep} lifetimes: max |MC sojourn - V^-1[:,s0]| = {zmax:.2e} "
    f"(entries ~ {Vinv_mc[:, s0].mean():.3f}); z-scores (batch-means SE): max |z| = {np.abs(zsc).max():.2f}, "
    f"rms z = {np.sqrt((zsc**2).mean()):.2f}; sum MC = {soj.sum():.4f} vs 1/gamma = {1 / gmc:.4f}")
RES["T1e_mc"] = dict(max_abs_err=float(zmax), mean_entry=float(Vinv_mc[:, s0].mean()), nrep=nrep,
                    max_abs_z=float(np.abs(zsc).max()), rms_z=float(np.sqrt((zsc**2).mean())))

log("-- Theorem 1'(f): Rayleigh/Ritz (symmetric generators)")
worst_ray, worst_modal, mono_fail, worst_full = 0.0, -1e9, 0, 0.0
for it in range(150):
    ins = random_instance(rng, "symmetric")
    M, q, b, g, n = ins["M"], ins["q"], ins["beta"], ins["gamma"], ins["n"]
    R0 = ng.R0_dense(M, q, b, g)
    # generalised symmetric eigenproblem
    A = np.diag(b * q)
    B = g * np.eye(n) - M.toarray()
    Rgen = sla.eigh(A, B, eigvals_only=True)[-1]
    worst_ray = max(worst_ray, abs(Rgen - R0) / R0)
    prev = -np.inf
    for nm in range(1, n + 1):
        Rr, diag_terms, mu, phi = ng.ritz_R0(M, q, b, g, nm)
        if Rr < prev - 1e-10 * R0:
            mono_fail += 1
        prev = Rr
    worst_full = max(worst_full, abs(Rr - R0) / R0)
    worst_modal = max(worst_modal, (diag_terms.max() - R0) / R0)  # must be <= 0
log(f"max rel |R0_generalised - R0_NGM|      = {worst_ray:.2e}")
log(f"max rel (largest diag. modal - R0)    = {worst_modal:.2e}   (<= 0 : modal formula is a LOWER bound)")
log(f"Ritz not monotone in #modes            : {mono_fail} failures;  full-basis rel err {worst_full:.2e}")
RES["T1f"] = dict(max_rel_rayleigh_err=worst_ray, max_rel_modal_minus_R0=worst_modal,
                  ritz_monotone_failures=mono_fail, full_basis_rel_err=worst_full)

# reversible non-uniform pi (departure model): symmetrised Rayleigh quotient
worst = 0.0
for it in range(100):
    ins = random_instance(rng, "departure")
    M, q, b, g, n, pi = ins["M"], ins["q"], ins["beta"], ins["gamma"], ins["n"], ins["pi"]
    sq = np.sqrt(pi)
    S = (M.toarray() * sq[None, :]) / sq[:, None]
    asym = np.abs(S - S.T).max()
    Rgen = sla.eigh(np.diag(b * q), g * np.eye(n) - 0.5 * (S + S.T), eigvals_only=True)[-1]
    worst = max(worst, abs(Rgen - ng.R0_dense(M, q, b, g)) / Rgen, asym)
log(f"departure model: max(asymmetry of Pi^-1/2 M Pi^1/2, rel err Rayleigh vs NGM) = {worst:.2e}")
RES["T1f_departure"] = worst

# ------------------------------------------------------------------ Theorem 2
log("== Theorem 2: beta<q>_pi/gamma <= R0 <= beta q_max/gamma, equality iff q constant ==")
N2 = 2000
viol_lo, viol_hi, min_gap_lo, min_gap_hi = 0.0, 0.0, np.inf, np.inf
pos = {k: [] for k in ["symmetric", "departure", "nonrev", "dense"]}
eq_err = 0.0
for it in range(N2):
    ins = random_instance(rng)
    M, q, b, g, pi = ins["M"], ins["q"], ins["beta"], ins["gamma"], ins["pi"]
    R0 = ng.R0_dense(M, q, b, g)
    lo, hi = ng.bounds(q, b, g, pi)
    viol_lo = max(viol_lo, (lo - R0) / R0)
    viol_hi = max(viol_hi, (R0 - hi) / R0)
    if q.max() - q.min() > 1e-3 * q.max():
        min_gap_lo = min(min_gap_lo, (R0 - lo) / R0)
        min_gap_hi = min(min_gap_hi, (hi - R0) / R0)
        pos[ins["kind"]].append((R0 - lo) / (hi - lo))
    # equality case
    qc = np.full(ins["n"], q.mean())
    eq_err = max(eq_err, abs(ng.R0_dense(M, qc, b, g) - b * q.mean() / g) / (b * q.mean() / g))
    # growth-rate bounds
    lam, *_ = ng.growth_rate_dense(M, q, b, g)
    viol_lo = max(viol_lo, ((b * pi @ q - g) - lam) / abs(g))
    viol_hi = max(viol_hi, (lam - (b * q.max() - g)) / abs(g))
log(f"instances {N2}; max relative violation of lower bound {viol_lo:.2e}, of upper bound {viol_hi:.2e}")
log(f"non-constant q: min relative gap to lower {min_gap_lo:.2e}, to upper {min_gap_hi:.2e} (strict)")
log(f"constant q: max rel |R0 - beta q/gamma| = {eq_err:.2e}")
for k, v in pos.items():
    v = np.array(v)
    log(f"  {k:10s}: n={len(v)}, normalised position (R0-lo)/(hi-lo): min {v.min():.3f} median {np.median(v):.3f} max {v.max():.3f}")
RES["T2"] = dict(n_instances=N2, max_viol_lower=viol_lo, max_viol_upper=viol_hi,
                 min_gap_lower=min_gap_lo, min_gap_upper=min_gap_hi, const_q_rel_err=eq_err,
                 position={k: [float(np.min(v)), float(np.median(v)), float(np.max(v))] for k, v in pos.items()})
np.save(OUT / "T2_positions.npy", np.array([(i, x) for i, k in enumerate(pos) for x in pos[k]]))

# ------------------------------------------------------------------ Theorem 3
log("== Theorem 3: mobility limits, monotonicity, asymptotics ==")
N3 = 300
Dgrid = np.logspace(-4, 4, 41)
mono_viol, lim_lo, lim_hi = 0.0, 0.0, 0.0
strict_fail = 0
exp_large, exp_small = [], []
conv_viol = 0.0
for it in range(N3):
    ins = random_instance(rng)
    M0, q, b, g, pi = ins["M"], ins["q"], ins["beta"], ins["gamma"], ins["pi"]
    # normalise the mobility scale so that the spectral gap of M0 is comparable to gamma
    R = np.array([ng.R0_dense(M0 * d, q, b, g) for d in Dgrid])
    d = np.diff(R)
    mono_viol = max(mono_viol, d.max() / R.mean())
    if q.max() - q.min() > 1e-3 * q.max() and not np.all(d < 0):
        strict_fail += 1
    lo, hi = ng.bounds(q, b, g, pi)
    Rinf = ng.R0_dense(M0 * 1e9, q, b, g)
    Rzero = ng.R0_dense(M0 * 1e-9, q, b, g)
    lim_lo = max(lim_lo, abs(Rinf - lo) / lo)
    lim_hi = max(lim_hi, abs(Rzero - hi) / hi)
    # growth rate: convex and non-increasing in D (perspective argument)
    lam = np.array([ng.growth_rate_dense(M0 * dd, q, b, g)[0] for dd in np.linspace(0.05, 5, 25)])
    conv_viol = max(conv_viol, -np.diff(lam, 2).min() / (abs(lam).max() + 1e-12), np.diff(lam).max() / (abs(lam).max() + 1e-12))
    # asymptotic expansions: error ratios when D changes by a factor 10 (should be ~100)
    C, qbar = ng.large_D_coefficient(M0, q, pi)
    evM = np.sort(np.linalg.eigvals(M0.toarray()).real)
    gapM = -evM[-2]                      # spectral gap of the unit-mobility generator
    rmax = float(np.abs(M0.diagonal()).max())  # largest exit rate
    sL = (b * q.max() + g) / gapM        # D >> sL  <=> mixing much faster than reaction
    sS = g / rmax                        # D << sS  <=> hopping much slower than recovery
    e1 = [abs(ng.R0_dense(M0 * D, q, b, g) - (b * qbar / g + b * C / (D * qbar))) for D in (30 * sL, 300 * sL)]
    # small D : hotspot zone
    Z = np.nonzero(q >= q.max() * (1 - 1e-12))[0]
    muZ = ng.mu_zone(M0, Z)
    e0 = [abs(ng.R0_dense(M0 * D, q, b, g) - hi * (1 - D * muZ / g)) for D in (1e-2 * sS, 1e-3 * sS)]
    if e1[1] > 1e-10 * lo and C > 1e-10:
        exp_large.append(e1[0] / e1[1])
    if e0[1] > 1e-10 * hi:
        exp_small.append(e0[0] / e0[1])
    if C < -1e-12:
        log("  NEGATIVE C", C, ins["kind"])
exp_large, exp_small = np.array(exp_large), np.array(exp_small)
log(f"instances {N3}; max increase of R0 between successive D (relative) = {mono_viol:.2e} (<=0 up to rounding)")
log(f"strict-decrease failures (non-constant q): {strict_fail}")
log(f"D->inf: max rel |R0 - beta<q>_pi/gamma| = {lim_lo:.2e};  D->0: max rel |R0 - beta q_max/gamma| = {lim_hi:.2e}")
log(f"growth rate lambda1(D): max violation of convexity / monotonicity = {conv_viol:.2e}")
log(f"large-D expansion error ratio err(30 s_L)/err(300 s_L): median {np.median(exp_large):.1f} (expect 100), "
    f"[{np.percentile(exp_large, 5):.1f}, {np.percentile(exp_large, 95):.1f}] (n={len(exp_large)})")
log(f"small-D expansion error ratio err(0.01 s_S)/err(0.001 s_S): median {np.median(exp_small):.1f} (expect 100), "
    f"[{np.percentile(exp_small, 5):.1f}, {np.percentile(exp_small, 95):.1f}] (n={len(exp_small)})")
RES["T3"] = dict(n_instances=N3, max_rel_increase=mono_viol, strict_failures=strict_fail,
                 limit_inf_rel_err=lim_lo, limit_zero_rel_err=lim_hi, convex_monotone_violation=conv_viol,
                 largeD_ratio_median=float(np.median(exp_large)), smallD_ratio_median=float(np.median(exp_small)))

# edge-wise monotonicity (symmetric rates) and zone lower bound (all generators)
log("-- Theorem 3': edge-wise monotonicity (symmetric) and zone lower bound")
edge_viol, zone_viol, n_edge, n_zone = 0.0, 0.0, 0, 0
dep_up, dep_down = 0, 0
for it in range(300):
    ins = random_instance(rng, "symmetric")
    lat, q, b, g = ins["lat"], ins["q"], ins["beta"], ins["gamma"]
    w = rng.lognormal(0, 1, len(lat.edges))
    M = ng.generator_from_edge_rates(lat.n, lat.edges, w, w)
    R0 = ng.R0_dense(M, q, b, g)
    for _ in range(3):
        w2 = w.copy()
        k = rng.integers(len(w))
        w2[k] *= rng.uniform(1.0, 20.0)
        R2 = ng.R0_dense(ng.generator_from_edge_rates(lat.n, lat.edges, w2, w2), q, b, g)
        edge_viol = max(edge_viol, (R2 - R0) / R0)  # must be <= 0
        n_edge += 1
    # departure model: increasing D at one cell can move R0 either way
    D = rng.lognormal(0, 0.8, lat.n)
    Rd = ng.R0_dense(ng.generator_departure(lat, D), q, b, g)
    D2 = D.copy()
    D2[rng.integers(lat.n)] *= 3.0
    Rd2 = ng.R0_dense(ng.generator_departure(lat, D2), q, b, g)
    if Rd2 > Rd * (1 + 1e-10):
        dep_up += 1
    elif Rd2 < Rd * (1 - 1e-10):
        dep_down += 1
for it in range(300):
    ins = random_instance(rng)
    M, q, b, g, n = ins["M"], ins["q"], ins["beta"], ins["gamma"], ins["n"]
    R0 = ng.R0_dense(M, q, b, g)
    for _ in range(3):
        Z = np.nonzero(rng.random(n) < rng.uniform(0.1, 0.9))[0]
        if len(Z) == 0:
            continue
        lb = b * q[Z].min() / (g + ng.mu_zone(M, Z))
        zone_viol = max(zone_viol, (lb - R0) / R0)
        n_zone += 1
log(f"symmetric: raising one edge rate, max rel increase of R0 over {n_edge} trials = {edge_viol:.2e} (<= 0)")
log(f"departure model: raising D at one cell -> R0 up in {dep_up}, down in {dep_down} of 300 trials (no monotonicity)")
log(f"zone bound beta q_Z/(gamma+mu_Z) <= R0: max rel violation over {n_zone} random zones = {zone_viol:.2e}")
RES["T3p"] = dict(edge_trials=n_edge, max_rel_increase_edge=edge_viol, departure_up=dep_up, departure_down=dep_down,
                  zone_trials=n_zone, max_rel_violation_zone=zone_viol)

# ------------------------------------------------------------------ Theorem 4
log("== Theorem 4: geometry; leaky room ==")
worst_uni = 0.0
for it in range(500):
    ins = random_instance(rng)
    q0 = rng.uniform(0.2, 3)
    R0 = ng.R0_dense(ins["M"], np.full(ins["n"], q0), ins["beta"], ins["gamma"])
    worst_uni = max(worst_uni, abs(R0 - ins["beta"] * q0 / ins["gamma"]) / R0)
log(f"(a) uniform q, any geometry/generator: max rel |R0 - beta q/gamma| = {worst_uni:.2e} (500 instances)")
RES["T4a_uniform_q_rel_err"] = worst_uni
# explicit size scan, reflecting square rooms, beta=0.5, gamma=0.14, q=1, D=1
size_scan = {}
for L in (3, 5, 10, 13, 20, 40):
    lat = ng.Lattice(np.ones((L, L), bool))
    M = ng.generator_symmetric(lat, np.ones(lat.n))
    R0 = ng.R0_sparse(M, np.ones(lat.n), 0.5, 0.14)
    lam = ng.growth_rate_sparse(M, np.ones(lat.n), 0.5, 0.14)
    size_scan[L] = dict(R0=R0, lambda1=lam)
    log(f"   L={L:3d}: R0 = {R0:.6f} (beta q/gamma = {0.5 / 0.14:.6f}), lambda1 = {lam:.6f}")
RES["T4a_size_scan"] = size_scan

worst_leak, worst_mu_ub, mono_mu, worst_het_lo, worst_het_hi = 0.0, 0.0, 0.0, 0.0, 0.0
worst_mu_lim = 0.0
for it in range(400):
    ins = random_instance(rng)
    M, q, b, g, n, pi = ins["M"], ins["q"], ins["beta"], ins["gamma"], ins["n"], ins["pi"]
    doors = np.nonzero(rng.random(n) < 0.25)[0]
    if len(doors) == 0 or len(doors) == n:
        doors = np.array([rng.integers(n)])
    kap = np.zeros(n)
    kap[doors] = rng.lognormal(0, 1.5, len(doors))
    mu = ng.mu_leak(M, kap)
    q0 = rng.uniform(0.2, 3)
    R0u = ng.R0_dense(M, np.full(n, q0), b, g, kappa=kap)
    worst_leak = max(worst_leak, abs(R0u - b * q0 / (g + mu)) / R0u)
    worst_mu_ub = max(worst_mu_ub, (mu - pi @ kap) / (pi @ kap))  # mu <= <kappa>_pi
    mu2 = ng.mu_leak(M, kap * 3)
    mono_mu = max(mono_mu, (mu - mu2) / mu)  # mu increasing in kappa
    # absorbing limit
    rest = np.setdiff1d(np.arange(n), doors)
    if len(rest):
        mu_abs = ng.mu_zone(M, rest)
        kbig = np.zeros(n)
        kbig[doors] = 1e9
        worst_mu_lim = max(worst_mu_lim, abs(ng.mu_leak(M, kbig) - mu_abs) / mu_abs, (mu - mu_abs) / mu_abs if mu > mu_abs else 0.0)
    # heterogeneous q with leak: nu = l*u normalised, quasi-stationary vectors of M - diag(kappa)
    s, u, l, _ = ng.perron_metzler(M.toarray() - np.diag(kap))
    nu = u * l / (u @ l)
    R0h = ng.R0_dense(M, q, b, g, kappa=kap)
    worst_het_lo = max(worst_het_lo, (b * (nu @ q) / (g + mu) - R0h) / R0h)
    worst_het_hi = max(worst_het_hi, (R0h - b * q.max() / (g + mu)) / R0h)
log(f"(b) leaky, uniform q: max rel |R0 - beta q/(gamma+mu_kappa)| = {worst_leak:.2e}")
log(f"    mu_kappa <= <kappa>_pi: max rel violation {worst_mu_ub:.2e};  mu increasing in kappa: max rel violation {mono_mu:.2e}")
log(f"    absorbing limit kappa->inf: max rel |mu - mu_Dirichlet| (and mu_kappa <= mu_Dirichlet) = {worst_mu_lim:.2e}")
log(f"(c) leaky, heterogeneous q: beta<q>_nu/(gamma+mu) <= R0 <= beta q_max/(gamma+mu): violations {worst_het_lo:.2e}, {worst_het_hi:.2e}")
RES["T4b"] = dict(uniform_q_rel_err=worst_leak, mu_upper_violation=worst_mu_ub, mu_monotone_violation=mono_mu,
                  absorbing_limit_err=worst_mu_lim, het_lower_violation=worst_het_lo, het_upper_violation=worst_het_hi)

# (d) 1-D room with absorbing ends: closed form mu = 2w(1-cos(pi/(L+1))) ~ pi^2 D/(L+1)^2
log("(d) 1-D corridor, both end cells open to the outside with exit rate w (= absorbing neighbours)")
one_d = {}
for L in (3, 5, 10, 13, 20, 40, 100):
    lat = ng.Lattice(np.ones((1, L), bool))
    w = 1.0
    M = ng.generator_symmetric(lat, np.full(L, w))
    kap = np.zeros(L)
    kap[0] = kap[-1] = w  # a jump out of either end leaves the room
    mu = ng.mu_leak(M, kap)
    exact = 2 * w * (1 - np.cos(np.pi / (L + 1)))
    cont = np.pi**2 * w / L**2
    R0 = ng.R0_dense(M, np.ones(L), 0.5, 0.14, kappa=kap)
    one_d[L] = dict(mu=mu, lattice_exact=exact, continuum_mu=cont, R0=R0, R0_with_continuum_mu=0.5 / (0.14 + cont))
    log(f"   L={L:3d}: mu = {mu:.6f}, 2w(1-cos(pi/(L+1))) = {exact:.6f}, pi^2 D/L^2 = {cont:.6f};"
        f" R0_room = {R0:.4f}, with pi^2 D/L^2: {0.5 / (0.14 + cont):.4f}")
RES["T4d_1d"] = one_d

# ------------------------------------------------------------------ Theorem 5
log("== Theorem 5: hotspot maps ==")
worst_early, worst_gen, worst_fd, worst_sum, worst_sym = 0.0, 0.0, 0.0, 0.0, 0.0
for it in range(200):
    ins = random_instance(rng)
    M, q, b, g, n = ins["M"], ins["q"], ins["beta"], ins["gamma"], ins["n"]
    q = np.maximum(q, 0.05)  # keep K primitive and q_r > 0 for the elasticity identity
    R0, w, z, K = ng.R0_dense(M, q, b, g, vectors=True)
    lam, u, l, _ = ng.growth_rate_dense(M, q, b, g)
    J = M.toarray() + np.diag(b * q - g)
    ev = np.linalg.eigvals(J)
    gap = lam - np.sort(ev.real)[-2]
    # (i) early-time profile from a point source
    T = 40.0 / max(gap, 1e-6)
    I0 = np.zeros(n)
    I0[rng.integers(n)] = 1.0
    IT = sla.expm((J - lam * np.eye(n)) * T) @ I0
    worst_early = max(worst_early, np.abs(IT / IT.sum() - u).max() / u.max())
    # (ii) generation profile
    x = I0.copy()
    for _ in range(3000):
        x = K @ x
        x /= x.sum()
    worst_gen = max(worst_gen, np.abs(x - w).max() / w.max())
    # (iii) elasticity
    e = ng.elasticity(w, z)
    worst_sum = max(worst_sum, abs(e.sum() - 1))
    r = rng.integers(n)
    h = 1e-5 * q[r]
    qp, qm = q.copy(), q.copy()
    qp[r] += h
    qm[r] -= h
    fd = (ng.R0_dense(M, qp, b, g) - ng.R0_dense(M, qm, b, g)) / (2 * h) * q[r] / R0
    worst_fd = max(worst_fd, abs(fd - e[r]) / max(e[r], 1e-12) if e[r] > 1e-6 else abs(fd - e[r]))
    if ins["kind"] == "symmetric":
        # e_r proportional to q_r x_r^2 with x the generalised eigenvector beta Q x = R0 V x
        vals, vecs = sla.eigh(np.diag(b * q), g * np.eye(n) - M.toarray())
        x = vecs[:, -1]
        e2 = q * x**2
        e2 /= e2.sum()
        worst_sym = max(worst_sym, np.abs(e2 - e).max() / e.max())
log(f"(i)  early-time profile -> Perron vector u of M+bQ-g: max rel err {worst_early:.2e}")
log(f"(ii) generation profile K^g -> w: max rel err {worst_gen:.2e}")
log(f"(iii) elasticities: |sum e - 1| <= {worst_sum:.2e}; finite-difference check max rel err {worst_fd:.2e}")
log(f"(iv) symmetric case e_r ~ q_r x_r^2: max rel err {worst_sym:.2e}")
RES["T5"] = dict(early_profile_err=worst_early, generation_profile_err=worst_gen, elasticity_sum_err=worst_sum,
                 elasticity_fd_err=worst_fd, symmetric_qx2_err=worst_sym)

# limits of the hotspot map
ins = random_instance(rng, "symmetric")
M0, q, b, g, n, pi = ins["M"], np.maximum(ins["q"], 0.05), ins["beta"], ins["gamma"], ins["n"], ins["pi"]
_, u_fast, _, _ = ng.growth_rate_dense(M0 * 1e7, q, b, g)
_, u_slow, _, _ = ng.growth_rate_dense(M0 * 1e-7, q, b, g)
Zmax = np.nonzero(q >= q.max() * (1 - 1e-12))[0]
log(f"(v) D->inf: max|u - pi| = {np.abs(u_fast - pi).max():.2e};  D->0: mass of u on argmax q = {u_slow[Zmax].sum():.6f}")
RES["T5_limits"] = dict(fast_err=float(np.abs(u_fast - pi).max()), slow_mass_on_argmax=float(u_slow[Zmax].sum()))

# ------------------------------------------------------------------ generalisations (remarks)
log("== Remarks: kernels / encounter-weighted q ==")
worst_kern_uni, viol_kmin, viol_kmax = 0.0, 0.0, 0.0
for it in range(300):
    ins = random_instance(rng)
    M, q, b, g, n = ins["M"], ins["q"], ins["beta"], ins["gamma"], ins["n"]
    C = rng.random((n, n)) * (rng.random((n, n)) < 0.4) + np.eye(n)
    C = C / C.sum(axis=0, keepdims=True)  # column-stochastic contact kernel
    Vinv = np.linalg.inv(g * np.eye(n) - M.toarray())
    for F in (b * C * q[None, :], b * q[:, None] * C):  # q at infector cell / at infectee cell
        R = np.abs(np.linalg.eigvals(F @ Vinv)).max()
        viol_kmin = max(viol_kmin, (b * q.min() / g - R) / max(R, 1e-12))
        viol_kmax = max(viol_kmax, (R - b * q.max() / g) / max(R, 1e-12))
    q0 = 1.3
    R = np.abs(np.linalg.eigvals(b * q0 * C @ Vinv)).max()
    worst_kern_uni = max(worst_kern_uni, abs(R - b * q0 / g) / R)
log(f"kernel models: uniform q -> R0 = beta q/gamma, max rel err {worst_kern_uni:.2e}; "
    f"beta q_min/gamma <= R0 <= beta q_max/gamma violations {viol_kmin:.2e}, {viol_kmax:.2e}")
RES["remarks_kernel"] = dict(uniform_rel_err=worst_kern_uni, viol_min=viol_kmin, viol_max=viol_kmax)

# ------------------------------------------------------------------ closed forms quoted in the text
log("== Closed forms quoted in the text ==")
latz = ng.Lattice(np.ones((15, 18), bool))
wz = 0.7
Mz = ng.generator_symmetric(latz, np.full(latz.n, wz))
zone_cf = []
for (y0, y1, x0, x1) in ((3, 8, 4, 11), (5, 7, 2, 12), (1, 13, 6, 9)):
    Zc = np.array([latz.index[y, x] for y in range(y0, y1) for x in range(x0, x1)])
    ly, lx = y1 - y0, x1 - x0
    muz = ng.mu_zone(Mz, Zc)
    exact = 2 * wz * (2 - np.cos(np.pi / (lx + 1)) - np.cos(np.pi / (ly + 1)))
    approx = wz * np.pi**2 * (1 / (lx + 1) ** 2 + 1 / (ly + 1) ** 2)
    zone_cf.append(dict(lx=lx, ly=ly, mu=muz, closed_form=exact, large_block_approx=approx))
    log(f"  interior block {lx}x{ly}: mu_Z = {muz:.12f}, closed form {exact:.12f}, pi^2 w[(lx+1)^-2+(ly+1)^-2] = {approx:.4f}")
Zc = np.array([latz.index[y, x] for y in range(0, 5) for x in range(4, 11)])  # block touching a wall in y
muz = ng.mu_zone(Mz, Zc)
exact = 2 * wz * (1 - np.cos(np.pi / (2 * 5 + 1))) + 2 * wz * (1 - np.cos(np.pi / 8))
log(f"  block 7x5 with one side on a wall: mu_Z = {muz:.12f}, closed form {exact:.12f}")
zone_cf.append(dict(wall_block_mu=muz, wall_block_closed_form=exact))
RES["zone_closed_forms"] = zone_cf
# symmetric rates: C = (1/n) dq^T (-M0)^+ dq  (Prop. 3.3(i))
worstC = 0.0
for _ in range(50):
    insC = random_instance(rng, "symmetric")
    M0c, qc_ = insC["M"].toarray(), insC["q"]
    Cgen, _ = ng.large_D_coefficient(M0c, qc_, insC["pi"])
    dq = qc_ - qc_.mean()
    Csym = dq @ np.linalg.pinv(-M0c) @ dq / insC["n"]
    worstC = max(worstC, abs(Cgen - Csym) / max(abs(Csym), 1e-12))
log(f"  symmetric rates: C = q^T x vs (1/n) dq^T (-M0)^+ dq, max rel. difference over 50 instances = {worstC:.2e}")
RES["C_symmetric_closed_form_rel_err"] = worstC
cfac = lambda N, p: 1 - (1 - (1 - p) ** N) / (N * p)
enc = {}
for N, ncell, qv in ((100, 234, 1.4), (100, 260, 1.4), (100, 400, 1.0), (400, 260, 1.4), (1600, 260, 1.4)):
    cval = cfac(N, 1 / ncell)
    enc[f"N={N},cells={ncell},q={qv}"] = dict(c=cval, R0_annealed=0.5 * qv / 0.14 * cval)
    log(f"  encounter factor c(N={N}, cells={ncell}) = {cval:.4f}; annealed R0 with q={qv}: {0.5 * qv / 0.14 * cval:.3f}")
RES["encounter_factors"] = enc

RES["runtime_s"] = time.time() - t0
log(f"total runtime {RES['runtime_s']:.1f} s")


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


(OUT / "verify_theorems.json").write_text(json.dumps(_clean(RES), indent=1))
(OUT / "verify_theorems.log").write_text("\n".join(LOG) + "\n")
