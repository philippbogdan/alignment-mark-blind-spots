"""Nystrom discretisation on Gauss-Legendre panels.

* far interactions: plain panel quadrature;
* self panel: kernel split K = A log|x-y| + B with the log-product weights of `quad.Tables`
  (A from Bessel functions, B diagonal limits from `MullerKernel.diag_B`);
* adjacent panels (sharing an endpoint): density interpolated to an auxiliary rule graded
  dyadically toward the shared endpoint (9 levels x 16 points).
Node coordinates for self / adjacent blocks are computed relative to a nearby point of the same
segment (no cancellation for tiny panels).
"""
import math
import numpy as np


# ---------------------------------------------------------------------------------- geometry helpers
def seg_rel(seg, mode, base, sig):
    """position(sig) - position(base) along segment, both distances measured from the segment
    start (mode 'start') or end (mode 'end'); stable for nearby points."""
    sig = np.asarray(sig, float)
    if seg.kind == "line":
        f = (sig - base) if mode == "start" else -(sig - base)
        return f * seg.tau[0], f * seg.tau[1]
    if mode == "start":
        phib = seg.phi0 + seg.sgn * base / seg.R
        dphi = seg.sgn * (sig - base) / seg.R
    else:
        phib = seg.phi1 - seg.sgn * base / seg.R
        dphi = -seg.sgn * (sig - base) / seg.R
    return seg._rel(phib, dphi)


def seg_abs(seg, mode, sig):
    s = sig if mode == "start" else seg.length - np.asarray(sig, float)
    return seg.at(s)


class Panel:
    def __init__(self, seg, mode, s0, s1):
        self.seg, self.mode, self.s0, self.s1 = seg, mode, float(s0), float(s1)
        self.length = abs(self.s1 - self.s0)
        self.speed = self.length / 2

    def sig(self, t):
        return self.s0 + (np.asarray(t) + 1) / 2 * (self.s1 - self.s0)

    def nodes(self, t, base=None, shift=(0.0, 0.0)):
        """dict of node data at parameters t.  Relative coordinates rx, rz: relative to the point
        at distance `base` (same segment/mode) plus `shift`; if base is None, rx = x (absolute)."""
        sg = self.sig(t)
        x, z, tx, tz, kap = seg_abs(self.seg, self.mode, sg)
        if base is None:
            rx, rz = x + shift[0], z + shift[1]
        else:
            rx, rz = seg_rel(self.seg, self.mode, base, sg)
            rx, rz = rx + shift[0], rz + shift[1]
        return dict(rx=np.asarray(rx, float), rz=np.asarray(rz, float), x=np.asarray(x, float),
                    z=np.asarray(z, float), nx=-np.asarray(tz, float), nz=np.asarray(tx, float),
                    kappa=np.asarray(kap, float))


def cat_nodes(lst):
    return {k: np.concatenate([d[k] for d in lst]) for k in lst[0]}


# ---------------------------------------------------------------------------------- assembly
def self_block(kern, tab, pan, base_nodes):
    """2x2 list of (n x n) weighted blocks for the self interaction of panel `pan`.
    base_nodes: node dict in a frame where coordinates are small (local)."""
    T = base_nodes
    M = kern.matrix(T, T)
    dx = T["rx"][:, None] - T["rx"][None, :]
    dz = T["rz"][:, None] - T["rz"][None, :]
    nt = len(T["rx"])
    nx = np.broadcast_to(T["nx"][:, None], (nt, nt))
    nz = np.broadcast_to(T["nz"][:, None], (nt, nt))
    mx = np.broadcast_to(T["nx"][None, :], (nt, nt))
    mz = np.broadcast_to(T["nz"][None, :], (nt, nt))
    A = kern.log_coeffs(dx, dz, nx, nz, mx, mz)
    Bd = kern.diag_B(T["nx"], T["nz"], T["kappa"])
    sp_ = pan.speed
    out = []
    idx = np.arange(nt)
    for a in range(4):
        blk = sp_ * (A[a] * tab.WL + (M[a] - A[a] * tab.logdt) * tab.w[None, :])
        blk[idx, idx] += sp_ * (A[a][idx, idx] * math.log(sp_) + Bd[a]) * tab.w
        out.append(blk)
    return out


def adjacent_block(kern, tab, tpan, spn, toward_minus, tnodes, aux_nodes_fn):
    """Weighted blocks for targets on panel tpan (node dict tnodes) and sources on panel spn whose
    endpoint t=-1 (toward_minus) or t=+1 is shared with tpan."""
    if toward_minus:
        ta, wa, I = tab.aux_t_minus, tab.aux_w, tab.I_minus
    else:
        ta, wa, I = tab.aux_t_plus, tab.aux_w_plus, tab.I_plus
    S = aux_nodes_fn(ta)
    M = kern.matrix(tnodes, S)
    wq = wa * spn.speed
    return [(m * wq[None, :]) @ I for m in M]


def assemble(kern, tab, panels, frame, adjacency, skip_blocks=()):
    """Full matrix K~ (2N x 2N, ordering [u nodes, v nodes]) for a list of panels.

    frame: 'abs' (coarse periodic mesh; self/adjacent blocks recomputed in local frames) or
           ('corner', corner_abs_point, sides) for RCIP local meshes, where every panel's
           coordinates are relative to the corner (base distance 0 on its segment end).
    adjacency: list of (p, q, toward_minus) ordered pairs.
    skip_blocks: set of (p, q) panel pairs whose blocks are left as plain quadrature (they are
           zeroed afterwards anyway).
    """
    n = tab.n
    npan = len(panels)
    N = n * npan
    if frame == "abs":
        nodes = [p.nodes(tab.t) for p in panels]
    else:
        nodes = [p.nodes(tab.t, base=0.0) for p in panels]
    ND = cat_nodes(nodes)
    w = np.concatenate([tab.w * p.speed for p in panels])
    M = kern.matrix(ND, ND)
    K = np.empty((2 * N, 2 * N), complex)
    K[:N, :N] = M[0] * w[None, :]
    K[:N, N:] = M[1] * w[None, :]
    K[N:, :N] = M[2] * w[None, :]
    K[N:, N:] = M[3] * w[None, :]
    del M

    def put(p, q, blks):
        rp = slice(p * n, (p + 1) * n)
        rq = slice(q * n, (q + 1) * n)
        rp2 = slice(N + p * n, N + (p + 1) * n)
        rq2 = slice(N + q * n, N + (q + 1) * n)
        K[rp, rq] = blks[0]
        K[rp, rq2] = blks[1]
        K[rp2, rq] = blks[2]
        K[rp2, rq2] = blks[3]

    for p, pan in enumerate(panels):
        if (p, p) in skip_blocks:
            continue
        if frame == "abs":
            loc = pan.nodes(tab.t, base=(pan.s0 + pan.s1) / 2)
        else:
            loc = nodes[p]
        put(p, p, self_block(kern, tab, pan, loc))
    for (p, q, tm) in adjacency:
        if (p, q) in skip_blocks:
            continue
        tp, sq = panels[p], panels[q]
        if frame == "abs" and tp.seg is sq.seg:
            base = sq.s0 if tm else sq.s1
            tn = tp.nodes(tab.t, base=base)
            fn = lambda ta, sq=sq, base=base: sq.nodes(ta, base=base)
        elif frame == "abs":
            # adjacent across a breakpoint (directly refined meshes): frame at the shared point
            assert tp.mode == sq.mode == "start"
            if tm:
                tn = tp.nodes(tab.t, base=tp.seg.length)
                fn = lambda ta, sq=sq: sq.nodes(ta, base=0.0)
            else:
                tn = tp.nodes(tab.t, base=0.0)
                fn = lambda ta, sq=sq: sq.nodes(ta, base=sq.seg.length)
        else:
            tn = nodes[p]
            fn = lambda ta, sq=sq: sq.nodes(ta, base=0.0)
        put(p, q, adjacent_block(kern, tab, tp, sq, tm, tn, fn))
    return K, ND, w
