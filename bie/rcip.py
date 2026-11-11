"""Coarse mesh with corner zones, and recursively compressed inverse preconditioning (RCIP,
Helsing 2008-2018) for every breakpoint of the profile.

For a breakpoint between segments A (before) and B (after), the coarse mesh carries two panels
of equal length s on each side (Gamma*).  The compressed inverse R (on Gamma*) is built by

    R_0 = (I + K_a(s 2^-n))^{-1}
    R_i = P_W^T ( F{R_{i-1}^{-1}} + I_b^o + K_b^o(s 2^-(n-i)) )^{-1} P,   i = 1..n

with type-b meshes of 6 panels (sizes s_i, s_i/2, s_i/2 per side), P the panel-to-children
interpolation and P_W = W_b P W_a^{-1}.  The coarse system is then (I + K^o R) rho~ = g~ and the
weight-corrected density is rho^ = R rho~.
"""
import math
import numpy as np
from .nystrom import Panel, assemble


# ---------------------------------------------------------------------------------- coarse mesh
def _interior(length, a, b, hmax, growth):
    """Breakpoints in (0, length) of a graded panel sequence, sizes ~a at 0, ~b at length, <= hmax."""
    if length <= 1e-14 * max(a, b):
        return np.array([0.0, length]) if length > 0 else np.array([0.0])
    m = 4001
    x = np.linspace(0.0, length, m)
    hloc = np.minimum.reduce([np.full(m, hmax), a + growth * x, b + growth * (length - x)])
    f = 1.0 / hloc
    F = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(x))])
    npan = max(1, int(math.ceil(F[-1] - 1e-9)))
    targets = np.linspace(0, F[-1], npan + 1)
    br = np.interp(targets, F, x)
    br[0], br[-1] = 0.0, length
    return br


def coarse_mesh(profile, hmax, growth=0.5, corner_factor=1.0, direct_levels=0):
    """Panels of the coarse mesh (absolute frame) and per-breakpoint info."""
    segs = profile.segs
    ns = len(segs)
    s_c = []
    for j in range(ns):
        a, b = segs[j - 1], segs[j]
        s_c.append(min(corner_factor * hmax, a.length / 4, b.length / 4))
    panels = []
    seg_panels = []
    for i, sg in enumerate(segs):
        a, b = s_c[i], s_c[(i + 1) % ns]
        L = sg.length
        inner = _interior(L - 2 * a - 2 * b, a, b, hmax, growth) + 2 * a
        br = np.concatenate([[0.0, a], inner, [L - b, L]])
        if direct_levels:
            # explicit dyadic refinement toward both breakpoints (validation of RCIP only)
            br = np.concatenate([br, a * 2.0 ** -np.arange(1, direct_levels + 1),
                                 L - b * 2.0 ** -np.arange(1, direct_levels + 1)])
        br = np.unique(br)
        ids = []
        for s0, s1 in zip(br[:-1], br[1:]):
            if s1 - s0 > 1e-12 * L:
                ids.append(len(panels))
                panels.append(Panel(sg, "start", s0, s1))
        seg_panels.append(ids)
    corners = []
    for j in range(ns):
        A = seg_panels[j - 1]
        B = seg_panels[j]
        gstar = [A[-2], A[-1], B[0], B[1]]
        if not direct_levels:
            for p in gstar:
                assert abs(panels[p].length - s_c[j]) < 1e-9 * s_c[j], (j, panels[p].length, s_c[j])
        corners.append(dict(j=j, segA=segs[j - 1], segB=segs[j], s=s_c[j], panels=gstar))
    return panels, corners


def cyclic_adjacency(panels, all_pairs=False):
    """Adjacent ordered pairs within the same segment (cross-segment pairs lie in Gamma*)."""
    adj = []
    n = len(panels)
    for p in range(n):
        q = (p + 1) % n
        if all_pairs or panels[p].seg is panels[q].seg:
            adj.append((p, q, True))
            adj.append((q, p, False))
    return adj


# ---------------------------------------------------------------------------------- RCIP
def _gstar(corner, sc):
    A, B = corner["segA"], corner["segB"]
    return [Panel(A, "end", 2 * sc, sc), Panel(A, "end", sc, 0.0),
            Panel(B, "start", 0.0, sc), Panel(B, "start", sc, 2 * sc)]


def _typeb(corner, sc):
    A, B = corner["segA"], corner["segB"]
    return [Panel(A, "end", 2 * sc, sc), Panel(A, "end", sc, sc / 2), Panel(A, "end", sc / 2, 0.0),
            Panel(B, "start", 0.0, sc / 2), Panel(B, "start", sc / 2, sc), Panel(B, "start", sc, 2 * sc)]


def _chain_adjacency(npan):
    adj = []
    for p in range(npan - 1):
        adj.append((p, p + 1, True))
        adj.append((p + 1, p, False))
    return adj


def _comp_idx(panel_ids, n, N):
    """indices of both unknown components for the given panels in a mesh of N nodes."""
    i = np.concatenate([np.arange(p * n, (p + 1) * n) for p in panel_ids])
    return np.concatenate([i, N + i])


def compute_R(kern, tab, corner, nlev):
    n = tab.n
    sc = corner["s"]
    # prolongation type-a (4 panels) -> type-b (6 panels), one component
    Pc = np.zeros((6 * n, 4 * n))
    Pc[0:n, 0:n] = np.eye(n)
    Pc[n:3 * n, n:2 * n] = tab.Pchild
    Pc[3 * n:5 * n, 2 * n:3 * n] = tab.Pchild
    Pc[5 * n:6 * n, 3 * n:4 * n] = np.eye(n)
    P = np.kron(np.eye(2), Pc)
    # weights (relative per level: all scale with sc, the ratio is level independent)
    wa = np.concatenate([tab.w * (sc / 2) * f for f in (1, 1, 1, 1)])
    wb = np.concatenate([tab.w * (sc / 2) * f for f in (1, 0.5, 0.5, 0.5, 0.5, 1)])
    wa2 = np.concatenate([wa, wa])
    wb2 = np.concatenate([wb, wb])
    PWT = (P * wb2[:, None]).T / wa2[:, None]          # W_a^{-1} P^T W_b
    Nb = 6 * n
    inner = _comp_idx([1, 2, 3, 4], n, Nb)
    adj_b = _chain_adjacency(6)
    adj_a = _chain_adjacency(4)
    skip = {(p, q) for p in (1, 2, 3, 4) for q in (1, 2, 3, 4)}
    # finest level
    s0 = sc * 2.0 ** (-nlev)
    Ka, _, _ = assemble(kern, tab, _gstar(corner, s0), "corner", adj_a)
    R = np.linalg.inv(np.eye(8 * n) + Ka)
    Io = np.eye(2 * Nb)
    Io[inner, inner] = 0.0
    for i in range(1, nlev + 1):
        si = sc * 2.0 ** (-(nlev - i))
        Kb, _, _ = assemble(kern, tab, _typeb(corner, si), "corner", adj_b, skip_blocks=skip)
        Kb[np.ix_(inner, inner)] = 0.0
        Mt = Io + Kb
        Mt[np.ix_(inner, inner)] += np.linalg.inv(R)
        R = PWT @ np.linalg.solve(Mt, P)
    return R
