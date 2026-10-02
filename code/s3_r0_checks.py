"""s3_r0_checks.py -- separate recomputation of every number quoted in Section 3
("The basic reproduction number of the mean-field model") that is not read directly from
data/plan_theory_checks.json, and a cross-check of those that are.

It uses only the scene definitions of the simulation part (code/bsc_sim/scenes/*.json through
bsc_sim.scenes / bsc_sim.theory.generator) and plain numpy; the formulas are those PROVED in
Section 3, in the article's index convention (column = source cell, M[y, x] = hop rate x -> y):

  * R0 = rho(beta Q (gamma 1 - D0 M0)^-1)                              (Theorem "ngm")
  * Rayleigh / modal formula and its truncations                      (Theorem "rayleigh", Corollary "modal")
  * fast-mixing coefficient  beta C / <q>,  C = q^T Y (Q pi - <q> pi),  Y = group inverse of -M0
    slow-mixing coefficient  mu_Z^0 = -s(M0[Z, Z]),  Z = argmax q      (Proposition "expansions")
  * zone bound beta q_Z / (gamma + mu_Z)                               (Proposition "zone")
  * offspring map n(x) = column sums of K                              (eq. "offspring")
  * two-zone room (q_W = 1.5, q_C = 0.5, <q> = 1.2)
  * three-cell contact-kernel counter-example, closed form (gamma + 4 D0)/(2 gamma + 6 D0)
                                                                       (Proposition "kernel")
  * R0 of the four canonical scenes against D0 (table of Section 3)
  * resolvent identity and both expansions on 200 random non-reversible generators

Output: ../data/s3_r0_checks.json      Run time: a few seconds, < 200 MB.
Run:    python s3_r0_checks.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                       # repository root
sys.path.insert(0, str(ROOT / "code/bsc_sim/src"))

from bsc_sim import scenes as sc             # noqa: E402
from bsc_sim import theory as th             # noqa: E402
from bsc_sim.scenes import Scene             # noqa: E402

BETA, GAMMA = 0.5, 0.14
OUT = HERE.parent / "data" / "s3_r0_checks.json"
PLAN = json.loads((HERE.parent / "data" / "plan_theory_checks.json").read_text())


def R0_of(M, q, P=None, gamma=GAMMA):
    """rho(F V^-1), F = beta P^T Q (beta Q if P is None), V = gamma 1 - M."""
    n = M.shape[0]
    K = (BETA * q)[:, None] * np.linalg.inv(gamma * np.eye(n) - M)
    if P is not None:
        K = P.T @ K
    return float(np.max(np.abs(np.linalg.eigvals(K)))), K


def spectral_abscissa(A):
    return float(np.max(np.linalg.eigvals(A).real))


def office():
    o = {}
    s = sc.load_scene("office", D0=1.0)
    n, q, zones = s.M, s.q, s.site_zone
    M0 = th.generator(s).toarray()                      # generator at D0 = 1 (symmetric rule)
    assert np.allclose(M0, M0.T) and np.allclose(M0.sum(axis=0), 0.0)
    pi = np.full(n, 1.0 / n)
    qbar, qmax = float(pi @ q), float(q.max())
    lower, upper = BETA * qbar / GAMMA, BETA * qmax / GAMMA
    o.update(cells=int(n), q_mean=qbar, q_max=qmax, q_min=float(q.min()), lower=lower, upper=upper)

    # --- R0 by the next-generation matrix and by the Rayleigh formula (symmetric eigenproblem)
    R0_1, K = R0_of(M0, q)
    B = GAMMA * np.eye(n) - M0                          # symmetric positive definite
    w, U = np.linalg.eigh(B)
    Bmh = (U / np.sqrt(w)) @ U.T                        # B^(-1/2)
    R0_rayleigh = float(np.linalg.eigvalsh(Bmh @ np.diag(BETA * q) @ Bmh)[-1])
    lam1 = spectral_abscissa(M0 + np.diag(BETA * q) - GAMMA * np.eye(n))
    lam1_rayleigh = float(np.linalg.eigvalsh(M0 + np.diag(BETA * q))[-1] - GAMMA)
    # unique root of s(M + beta Q / R - gamma 1) = 0  (Theorem ngm (b))
    root_residual = spectral_abscissa(M0 + np.diag(BETA * q) / R0_1 - GAMMA * np.eye(n))
    o["D0=1"] = dict(R0=R0_1, R0_rayleigh=R0_rayleigh, lambda1=lam1, lambda1_rayleigh=lam1_rayleigh,
                     root_residual=root_residual, doubling_time=float(np.log(2) / lam1))

    # --- offspring map: column sums of K
    nx_ = K.sum(axis=0)
    o["offspring"] = dict(min=float(nx_.min()), max=float(nx_.max()), pi_mean=float(pi @ nx_),
                          beta_qmin_over_gamma=BETA * float(q.min()) / GAMMA)

    # --- modal formula (Corollary modal): eigenpairs of -M0, coupling matrix Lambda
    mu, phi = np.linalg.eigh(-M0)
    mu[0] = 0.0
    Lam = phi.T @ (q[:, None] * phi)
    g = 1.0 / np.sqrt(GAMMA + mu)
    Rm = BETA * Lam * g[:, None] * g[None, :]
    diag = np.diag(Rm)
    o["modal"] = dict(diagonal_max=float(diag.max()), diagonal_argmax=int(np.argmax(diag)),
                      diagonal_0=float(diag[0]),
                      truncations={str(m): float(np.linalg.eigvalsh(Rm[:m, :m])[-1]) for m in (1, 2, 3, 5, 10, 20, 50, n)},
                      mu_1=float(mu[1]))

    # --- expansions (Proposition expansions)
    Pproj = np.outer(pi, np.ones(n))
    Y = np.linalg.inv(Pproj - M0) - Pproj               # group inverse of -M0
    y = Y @ (q * pi - qbar * pi)
    assert abs(y.sum()) < 1e-10 and np.allclose(-M0 @ y, q * pi - qbar * pi)
    C = float(q @ y)
    dq = q - qbar
    C_sym = float(dq @ np.linalg.pinv(-M0) @ dq) / n    # symmetric closed form
    Z = np.nonzero(np.isclose(q, qmax))[0]
    muZ0 = -spectral_abscissa(M0[np.ix_(Z, Z)])
    o["expansions"] = dict(C=C, C_symmetric_form=C_sym, fast_coefficient=BETA * C / qbar,
                           muZ0=muZ0, slow_slope=muZ0 / GAMMA, Z_cells=int(len(Z)))
    rows = []
    for D0 in (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4):
        R = R0_of(D0 * M0, q)[0]
        fast = lower + BETA * C / (qbar * D0)
        slow = upper * (1 - D0 * muZ0 / GAMMA)
        rows.append(dict(D0=D0, R0=R, fast=fast, slow=slow, err_fast=R - fast, err_slow=R - slow))
    o["expansions"]["rows"] = rows

    # --- zone bound, Z = meeting room (zone code M), D0 = 1
    ZM = np.nonzero(zones == "M")[0]
    muM = -spectral_abscissa(M0[np.ix_(ZM, ZM)])
    o["zone_meeting_room"] = dict(cells=int(len(ZM)), q_Z=float(q[ZM].min()), mu_Z=muM,
                                  exit_time_days=1.0 / muM, bound=BETA * float(q[ZM].min()) / (GAMMA + muM))

    # --- mobility sweep (monotone decrease), including the physical values of Section 2
    sweep = []
    for D0 in (0.01, 0.1, 1.0, 10.0, 100.0, 2.4e3, 4.87e4):
        R = R0_of(D0 * M0, q)[0]
        sweep.append(dict(D0=D0, R0=R, excess_pct=100 * (R / lower - 1)))
    o["sweep"] = sweep
    return o


def two_zone():
    """Two-zone room, as built in
    code/bsc_sim/experiments/05_heterogeneity.py: 20 x 13 cells, central corridor band of 6 columns
    (30 % of the area) with D = 2 D0 and q_C = 0.5, workspace with D = 0.3 D0 and q_W = 1.5."""
    nx, ny = 20, 13
    zone = np.full((ny, nx), "W")
    zone[:, 7:13] = "C"
    access = np.ones((ny, nx), dtype=bool)
    q = np.where(zone == "W", 1.5, 0.5)
    D = np.where(zone == "W", 0.3, 2.0)
    s = Scene(name="twozone", nx=nx, ny=ny, a_m=1.5, zone=zone, access=access, D_grid=D, q_grid=q, N=100, meta={})
    M0 = th.generator(s).toarray()
    qv = s.q
    out = dict(q_mean=float(qv.mean()), lower_over_beta_gamma=float(qv.mean()), upper_over_beta_gamma=float(qv.max()),
               rows=[])
    for D0 in (1.0, 10.0, 100.0, 1000.0):
        R = R0_of(D0 * M0, qv)[0]
        out["rows"].append(dict(D0=D0, R0=R, R0_over_beta_gamma=R * GAMMA / BETA))
    return out


def kernel_example():
    P = np.array([[0.5, 0.5, 0.0], [1 / 3, 1 / 3, 1 / 3], [0.0, 0.5, 0.5]])
    M0 = np.array([[-1.0, 1.0, 0.0], [1.0, -2.0, 1.0], [0.0, 1.0, -1.0]])
    out = {}
    for qv in ([1.0, 0.0, 1.0], [1.0, 0.2, 1.0]):
        qv = np.array(qv)
        rows = []
        for D0 in (0.0, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 1e3):
            R = R0_of(D0 * M0, qv, P=P)[0] * GAMMA / BETA
            row = dict(D0=D0, R0_over_beta_gamma=R)
            if qv[1] == 0.0:
                row["closed_form"] = (GAMMA + 4 * D0) / (2 * GAMMA + 6 * D0)
            rows.append(row)
        out[f"q={qv.tolist()}"] = dict(mean_q=float(qv.mean()), rows=rows,
                                       slow_limit_rho_PTQ=float(np.max(np.abs(np.linalg.eigvals(P.T @ np.diag(qv))))))
    # same-cell kernel on the same three cells: the bounds of Theorem "bounds" hold
    rows = []
    qv = np.array([1.0, 0.0, 1.0])
    for D0 in (0.0, 1e-2, 1.0, 1e3):
        rows.append(dict(D0=D0, R0_over_beta_gamma=R0_of(D0 * M0, qv)[0] * GAMMA / BETA))
    out["same_cell_q=[1,0,1]"] = rows
    return out


def scenes_table():
    """Table of Section 3: R0 of the four canonical scenes against the mobility scale, recomputed here
    (1 m contact kernel and same-cell kernel) and compared with data/plan_theory_checks.json."""
    out = {}
    worst = 0.0
    for name in ("office", "supermarket", "classroom", "metro"):
        s = sc.load_scene(name, D0=1.0)
        M0 = th.generator(s).toarray()
        P = th.contact_kernel(s, th.kernel_radius_cells(s, 1.0)).toarray()
        q = s.q
        lower, upper = BETA * float(q.mean()) / GAMMA, BETA * float(q.max()) / GAMMA
        row = dict(cells=int(s.M), kernel_cells_mean=float((P > 0).sum(axis=1).mean()), lower=lower, upper=upper,
                   slow_limit_1m=float(np.max(np.abs(np.linalg.eigvals(P.T * (BETA * q)[None, :])))) / GAMMA)
        for D0 in (1.0, 10.0, 100.0, 2400.0):
            R1m = R0_of(D0 * M0, q, P=P)[0]
            Rsc = R0_of(D0 * M0, q)[0]
            row[f"D0={D0:g}"] = dict(R0_1m=R1m, R0_samecell=Rsc, excess_1m_pct=100 * (R1m / lower - 1),
                                     excess_samecell_pct=100 * (Rsc / lower - 1))
            pl = PLAN["scenes"][name]
            worst = max(worst, abs(R1m - pl[f"R0_1m_D0={D0:g}"]), abs(Rsc - pl[f"R0_samecell_D0={D0:g}"]))
        out[name] = row
    out["max_abs_diff_to_plan_theory_checks"] = worst
    assert worst < 1e-8, worst
    return out


def nonreversible_check(n_inst=200, seed=20261001):
    """Proposition "expansions" and the resolvent identity in the article's convention (column = source),
    on random NON-reversible generators with several cells sharing the maximal q.
    For each instance the error of the fast (slow) expansion is evaluated at two mobilities a decade
    apart; the fast error must fall about 100-fold (O(D0^-2) remainder) and the slow error divided by D0
    must fall (o(D0) remainder)."""
    rng = np.random.default_rng(seed)
    worst_res, ratios_fast, ratios_slow, minC, n_floor = 0.0, [], [], np.inf, 0
    for _ in range(n_inst):
        n = int(rng.integers(4, 10))
        W = rng.lognormal(0, 1, (n, n)) * (rng.random((n, n)) < 0.6)
        np.fill_diagonal(W, 0.0)
        for i in range(n):
            W[(i + 1) % n, i] += rng.uniform(0.1, 1.0)      # a cycle makes the generator irreducible
        M0 = W - np.diag(W.sum(axis=0))                      # M0[y, x] = rate x -> y, zero column sums
        q = rng.uniform(0.2, 1.0, n)
        k = int(rng.integers(1, 4))
        Z = rng.choice(n, size=k, replace=False)
        q[Z] = 1.3                                           # |Z| = k cells share q_max
        ev, vec = np.linalg.eig(M0)
        pi = np.abs(vec[:, np.argmin(np.abs(ev))].real)
        pi /= pi.sum()
        qbar = float(pi @ q)
        E0 = np.outer(pi, np.ones(n))
        Y = np.linalg.inv(E0 - M0) - E0
        eps = 1e-3                                           # resolvent identity at D0 = 1000
        lhs = np.linalg.inv(GAMMA * np.eye(n) - M0 / eps)
        rhs = E0 / GAMMA + eps * Y @ np.linalg.inv(np.eye(n) + eps * GAMMA * Y)
        worst_res = max(worst_res, float(np.abs(lhs - rhs).max()))
        C = float(q @ Y @ (q * pi))
        minC = min(minC, C)
        muZ0 = -spectral_abscissa(M0[np.ix_(Z, Z)])
        lower, upper = BETA * qbar / GAMMA, BETA * 1.3 / GAMMA
        ef = [abs(R0_of(D * M0, q)[0] - lower - BETA * C / (qbar * D)) for D in (1e2, 1e3)]
        es = [abs(R0_of(D * M0, q)[0] - upper * (1 - D * muZ0 / GAMMA)) / D for D in (1e-3, 1e-4)]
        if ef[0] > 1e-8:                                     # well above the rounding floor (about 1e-11) of the eigen-solver
            ratios_fast.append(ef[0] / ef[1])
        else:
            n_floor += 1
        ratios_slow.append(es[0] / max(es[1], 1e-300))
    return dict(n_instances=n_inst, resolvent_identity_max_abs_err=worst_res, min_C=float(minC),
                fast_instances_at_rounding_floor_excluded=n_floor,
                fast_error_ratio_per_decade_median=float(np.median(ratios_fast)),
                fast_error_ratio_per_decade_min=float(np.min(ratios_fast)),
                slow_error_over_D0_ratio_per_decade_median=float(np.median(ratios_slow)),
                slow_error_over_D0_ratio_per_decade_min=float(np.min(ratios_slow)))


def main():
    res = dict(beta=BETA, gamma=GAMMA, beta_over_gamma=BETA / GAMMA, office=office(), two_zone=two_zone(),
               kernel_example=kernel_example(), scenes=scenes_table(), nonreversible_check=nonreversible_check())
    # cross-check against plan_theory_checks.json
    po = PLAN["office"]
    chk = {
        "R0_D0=1": (res["office"]["D0=1"]["R0"], po["D0=1"]["R0"]),
        "lambda1_D0=1": (res["office"]["D0=1"]["lambda1"], po["D0=1"]["lambda1"]),
        "fast_coefficient": (res["office"]["expansions"]["fast_coefficient"], po["largeD_coefficient_betaC_over_qbar"]),
        "slow_slope": (res["office"]["expansions"]["slow_slope"], po["smallD_slope_muZ0_over_gamma"]),
        "diagonal_modal_max": (res["office"]["modal"]["diagonal_max"], po["diagonal_modal_max"]),
        "galerkin_3": (res["office"]["modal"]["truncations"]["3"], po["galerkin_R0_by_modes"]["3"]),
        "galerkin_10": (res["office"]["modal"]["truncations"]["10"], po["galerkin_R0_by_modes"]["10"]),
        "zone_bound": (res["office"]["zone_meeting_room"]["bound"], po["zone_bound_meeting_room"]["bound"]),
        "offspring_min": (res["office"]["offspring"]["min"], po["D0=1"]["offspring_min"]),
        "offspring_max": (res["office"]["offspring"]["max"], po["D0=1"]["offspring_max"]),
    }
    res["crosscheck_with_plan_theory_checks"] = {k: dict(here=a, plan_theory_checks=b, abs_diff=abs(a - b)) for k, (a, b) in chk.items()}
    worst = max(v["abs_diff"] for v in res["crosscheck_with_plan_theory_checks"].values())
    res["crosscheck_max_abs_diff"] = worst
    assert worst < 1e-8, worst
    kc = res["kernel_example"]["q=[1.0, 0.0, 1.0]"]["rows"]
    res["kernel_closed_form_max_abs_diff"] = max(abs(r["R0_over_beta_gamma"] - r["closed_form"]) for r in kc)
    assert res["kernel_closed_form_max_abs_diff"] < 1e-12
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
