#!/usr/bin/env python3
"""NOVA-D population calibration and feasibility solver (deterministic; no selector, no simulation).

For each grid condition, solves the DGP structural parameters so that the POPULATION oracle
quantities hit their targets exactly, then verifies them against closed forms.

Family (NOVA_D_PROTOCOL_v0.2 section 4), with g = G_(1) = 1 by scale invariance:
    G_1 = 1,  G_2 = 1 - Gamma_rel,  G_j = (1 - Gamma_rel) * r^(j-2),  j = 3..p
    H on exactly one designated candidate; H_j = 0 otherwise
    R_empty = sum_j G_j + sum_j H_j + d_T * sigma_eps^2

Exact identities used (squared-error Bayes risk; Pythagoras gives R_S = R_empty - V_S):
    A        = V_{j_G} / R_empty                     (j_G = argmax G)
    rho_pop  = (V_{j*}/V_{j_G} - 1) / (1 - A)        (j* = argmax V)
Both are EXACT for this family, not approximations.

Unknowns per condition: r (G-profile decay), H, sigma_eps^2.
Equations: the A target, the rho_pop target, and a prespecified irreducible-noise share
    psi = d_T * sigma_eps^2 / R_empty.
"""
from __future__ import annotations
import json
from dataclasses import dataclass, asdict

P = 20          # candidates
D_T = 20        # targets
TOL = 1e-10


def sumG(r: float, gam: float, p: int = P) -> float:
    """sum_j G_j with G_1 = 1, G_2 = 1-gam, G_j = (1-gam) r^(j-2)."""
    n = p - 1                                    # terms j=2..p
    geo = n if abs(r - 1.0) < 1e-15 else (1.0 - r ** n) / (1.0 - r)
    return 1.0 + (1.0 - gam) * geo


def G_at(rank: int, r: float, gam: float) -> float:
    if rank == 1:
        return 1.0
    return (1.0 - gam) * r ** (rank - 2)


def required_H(mech_rank: int, r: float, gam: float, A: float, rho: float) -> float:
    """H on the designated candidate so that it is j* with the target rho_pop.

    j_G is rank 1 (G_1 = 1 is the max). V_{j*} = G_{rank} + H, V_{j_G} = 1.
    rho = (V_{j*}/1 - 1)/(1-A)  =>  G_rank + H - 1 = rho (1-A).
    """
    return 1.0 - G_at(mech_rank, r, gam) + rho * (1.0 - A)


def solve_r(mech_rank: int, gam: float, A: float, rho: float, psi: float,
            lo: float = 1e-9, hi: float = 0.999999) -> float | None:
    """Solve (1-psi)/A = sumG(r) + H(r) for r by bisection. None if infeasible."""
    def f(r):
        return sumG(r, gam) + required_H(mech_rank, r, gam, A, rho) - (1.0 - psi) / A
    flo, fhi = f(lo), f(hi)
    if flo > 0:            # even maximal concentration leaves too much value spread
        return None
    if fhi < 0:            # even maximal spread cannot absorb the budget
        return None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) <= 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass
class Cell:
    mech: str
    rank: int          # G-rank of the designated nonlinear candidate (1 = none/null)
    gam: float
    A_target: float
    rho_target: float
    psi: float
    feasible: bool
    r: float | None = None
    H: float | None = None
    sigma2: float | None = None
    A_oracle: float | None = None
    rho_oracle: float | None = None
    noise_share: float | None = None


def build(mech: str, rank: int, gam: float, A: float, rho: float, psi: float) -> Cell:
    c = Cell(mech, rank, gam, A, rho, psi, feasible=False)
    if rho == 0.0:                      # M1 / M3: H identically zero
        r = solve_r(1, gam, A, 0.0, psi)
        if r is None:
            return c
        H = 0.0
    else:
        r = solve_r(rank, gam, A, rho, psi)
        if r is None:
            return c
        H = required_H(rank, r, gam, A, rho)
        if H <= 0:
            return c
    SG = sumG(r, gam)
    R0 = 1.0 / A                        # since V_{j_G} = G_1 = 1 and A = 1/R0
    s2 = (R0 - SG - H) / D_T
    if s2 < -TOL:
        return c
    s2 = max(s2, 0.0)
    # closed-form verification of the oracle quantities
    R0_chk = SG + H + D_T * s2
    V_jG = 1.0
    V_star = max(G_at(k, r, gam) + (H if (rho > 0 and k == rank) else 0.0)
                 for k in range(1, P + 1))
    A_or = V_jG / R0_chk
    rho_or = (V_star / V_jG - 1.0) / (1.0 - A_or)
    if rho == 0.0 and abs(rho_or) > 1e-9:
        return c
    if rho > 0:
        # the designated candidate must BE the argmax of V, and sit at its stated G-rank
        v = [G_at(k, r, gam) + (H if k == rank else 0.0) for k in range(1, P + 1)]
        if v.index(max(v)) + 1 != rank:
            return c
        gs = sorted(range(1, P + 1), key=lambda k: -G_at(k, r, gam))
        if gs.index(rank) + 1 != rank:
            return c
    c.feasible = (abs(A_or - A) < 1e-9 and abs(rho_or - rho) < 1e-9)
    c.r, c.H, c.sigma2 = r, H, s2
    c.A_oracle, c.rho_oracle = A_or, rho_or
    c.noise_share = D_T * s2 / R0_chk
    return c


MECHS = [("M1", 1, 0.20, [0.0]),
         ("M3", 1, 0.02, [0.0]),
         ("M2", 2, 0.02, [0.05, 0.10, 0.40]),
         ("M4", 7, 0.20, [0.05, 0.10, 0.40])]


def sweep(A_grid, psi):
    rows = []
    for mech, rank, gam, rhos in MECHS:
        for rho in rhos:
            for A in A_grid:
                rows.append(build(mech, rank, gam, A, rho, psi))
    return rows


def ceiling(mech, rank, gam, rho, psi, step=0.0005):
    """Largest A on a fine grid for which the condition is feasible."""
    best = None
    A = step
    while A < 0.95:
        if build(mech, rank, gam, A, rho, psi).feasible:
            best = A
        A += step
    return best


if __name__ == "__main__":
    PSI = 0.25     # prespecified minimum irreducible-noise share of base-panel risk
    print(f"irreducible-noise share psi = {PSI}\n")

    print("=== owner's confirmed A grid {0.05, 0.15, 0.40} ===")
    rows = sweep([0.05, 0.15, 0.40], PSI)
    nf = [r for r in rows if not r.feasible]
    print(f"conditions: {len(rows)}   feasible: {len(rows)-len(nf)}   INFEASIBLE: {len(nf)}")
    for r in nf:
        print(f"  INFEASIBLE  {r.mech}  rho_pop={r.rho_target:<4}  A={r.A_target}")

    print("\n=== structural ceiling on A, per mechanism and advantage ===")
    print(f"{'mech':6}{'rho_pop':>9}{'max feasible A':>16}")
    for mech, rank, gam, rhos in MECHS:
        for rho in rhos:
            c = ceiling(mech, rank, gam, rho, PSI)
            print(f"{mech:6}{rho:>9}{(f'{c:.3f}' if c else 'none'):>16}")

    print("\n=== candidate replacement grid {0.05, 0.15, 0.30} ===")
    rows2 = sweep([0.05, 0.15, 0.30], PSI)
    nf2 = [r for r in rows2 if not r.feasible]
    print(f"conditions: {len(rows2)}   feasible: {len(rows2)-len(nf2)}   INFEASIBLE: {len(nf2)}")
    for r in nf2:
        print(f"  INFEASIBLE  {r.mech}  rho_pop={r.rho_target:<4}  A={r.A_target}")

    print("\n=== calibrated parameters, primary and control conditions ===")
    for mech, rank, gam, rho, A, lab in [("M2", 2, 0.02, 0.10, 0.15, "PRIMARY"),
                                         ("M2", 2, 0.02, 0.40, 0.15, "positive control"),
                                         ("M4", 7, 0.20, 0.10, 0.15, "aim control"),
                                         ("M3", 1, 0.02, 0.00, 0.15, "instability control")]:
        c = build(mech, rank, gam, A, rho, PSI)
        print(f"{lab:20} {mech}  A={A} rho={rho}  ->  r={c.r:.6f}  H={c.H:.6f}  "
              f"sigma_eps^2={c.sigma2:.6f}  noise={c.noise_share:.3f}  "
              f"[verify A={c.A_oracle:.10f} rho={c.rho_oracle:.10f}]")

    with open("nova_d_calibration_sweep.json", "w") as fh:
        json.dump({"psi": PSI,
                   "owner_grid": [asdict(r) for r in rows],
                   "replacement_grid": [asdict(r) for r in rows2]}, fh, indent=1)
    print("\nwrote nova_d_calibration_sweep.json")
