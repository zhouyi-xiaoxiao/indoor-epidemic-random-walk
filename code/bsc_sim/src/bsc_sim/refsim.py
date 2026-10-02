"""Readable pure-Python reference implementation of the simulator.

Same model and same exact (Gillespie direct-method) algorithm as the C core
in ``gillespie.c``; slow, but short enough to read in one sitting.  It is
used (i) to cross-validate the C core and (ii) as a readable listing of the
algorithm.

State
    pos[i]    cell of person i            (0 .. M-1, walkable cells only)
    state[i]  0 susceptible, 1 infective, 2 removed

Events and rates
    hop        person i: x -> neighbour y          W(x -> y)
    recovery   infective i                         gamma
    infection  susceptible j at y by infective i at x
               pair-hazard rule:  beta q(x) P(y|x) / rho_bar
               per-cell rule   :  beta q(y) / n_total(y)   if x == y, else 0
"""
from __future__ import annotations

import numpy as np

from .scenes import Scene
from .theory import contact_kernel, hop_rates


def simulate_one(scene: Scene, beta: float, gamma: float, rng,
                 kernel_radius: float = 0.0, mode: str = "kernel",
                 N: int | None = None, rho_ref: float | None = None,
                 index_site: int | None = None, gmax: int = -1,
                 t_max: float = 150.0):
    """Run one epidemic; return dict(final_size, index_offspring, t_inf, ...)."""
    M = scene.M
    N = scene.N if N is None else N
    ptr, idx, W = hop_rates(scene)
    Wtot = np.add.reduceat(W, ptr[:-1])
    P = contact_kernel(scene, kernel_radius).toarray()
    rho = (N - 1) / M if rho_ref is None else rho_ref
    q = scene.q

    pos = rng.integers(0, M, size=N)
    if index_site is not None:
        pos[0] = index_site
    state = np.zeros(N, dtype=int)
    state[0] = 1
    gen = -np.ones(N, dtype=int)
    gen[0] = 0
    infector = -np.ones(N, dtype=int)
    t_inf = -np.ones(N)
    t_inf[0] = 0.0
    t = 0.0
    n_events = 0

    def pair_hazard(i, j):
        """Rate at which infective i infects susceptible j."""
        x, y = pos[i], pos[j]
        if mode == "percell":
            return beta * q[y] / np.sum(pos == y) if x == y else 0.0
        return beta * q[x] * P[x, y] / rho

    while True:
        infectives = np.flatnonzero(state == 1)
        if len(infectives) == 0 or t >= t_max:
            break
        susceptibles = np.flatnonzero(state == 0)
        hop = Wtot[pos]                                   # (N,)
        rec = gamma * np.ones(len(infectives))
        inf = np.array([[pair_hazard(i, j) for j in susceptibles]
                        for i in infectives]).reshape(len(infectives),
                                                      len(susceptibles))
        total = hop.sum() + rec.sum() + inf.sum()
        t += rng.exponential(1.0 / total)
        if t >= t_max:
            break
        u = rng.random() * total
        n_events += 1
        if u < hop.sum():                                 # ---- hop
            i = np.searchsorted(np.cumsum(hop), u, side="right")
            x = pos[i]
            w = W[ptr[x]:ptr[x + 1]]
            k = rng.choice(len(w), p=w / w.sum())
            pos[i] = idx[ptr[x] + k]
        elif u < hop.sum() + rec.sum():                   # ---- recovery
            i = infectives[rng.integers(len(infectives))]
            state[i] = 2
        else:                                             # ---- infection
            k = rng.choice(inf.size, p=inf.ravel() / inf.sum())
            i = infectives[k // len(susceptibles)]
            j = susceptibles[k % len(susceptibles)]
            gen[j] = gen[i] + 1
            infector[j] = i
            t_inf[j] = t
            state[j] = 2 if (gmax >= 0 and gen[j] > gmax) else 1

    return dict(final_size=int((t_inf >= 0).sum()),
                index_offspring=int((infector == 0).sum()),
                t_inf=t_inf, gen=gen, infector=infector, t_end=t,
                n_events=n_events)
