"""Finite-element solver (NGSolve) for the alignment mark of MODEL.md.

Shares no code with bie/: its own profile construction (straight edges and circular arcs, curved
elements), periodic mesh in x, geometric (hp) refinement at the silicon corners, and exact
Dirichlet-to-Neumann (Fourier-mode) radiation conditions.  solve_rg computes the reflection
matrix of the grating (silicon below the profile, an oxide half-space above) at a plane above
the line; cover adds the oxide/air interface by the exact multiple-reflection formula, so the
oxide thickness enters analytically.  readings returns what a channel measures.

TE (u = E_y):  -div grad u - k0^2 eps u = 0
TM (u = H_y):  -div (1/eps) grad u - k0^2 u = 0
"""
import math
import numpy as np

from . import geometry as G
from .materials import eps_si, eps_sio2

_TM_ACTIVE = [False]


def _threaded(fn):
    """Run fn inside NGSolve's TaskManager (multi-threaded assembly and factorisation; about 2x
    faster on 4 cores, identical results), unless one is already active."""
    import functools

    @functools.wraps(fn)
    def wrapper(*args, **kw):
        if _TM_ACTIVE[0]:
            return fn(*args, **kw)
        import os
        import ngsolve as ng
        nthr = int(os.environ.get("OMP_NUM_THREADS") or 0)
        if nthr > 0:
            ng.SetNumThreads(nthr)
        _TM_ACTIVE[0] = True
        try:
            with ng.TaskManager():
                return fn(*args, **kw)
        finally:
            _TM_ACTIVE[0] = False
    return wrapper

P = G.P


class _Solver:
    """Sparse factorisation of the assembled NGSolve matrix restricted to free dofs, with one
    step of iterative refinement per solve (residual computed with the exact sparse matrix)."""

    def __init__(self, amat, fes, idx, kind="sparsecholesky"):
        import ngsolve as ng
        self.amat, self.fes, self.idx = amat, fes, idx
        self.inv = amat.Inverse(fes.FreeDofs(), inverse=kind)
        self.gv = ng.GridFunction(fes).vec
        self.rv = self.gv.CreateVector()
        self.sv = self.gv.CreateVector()
        self.maxres = 0.0

    def _apply(self, b):
        self.rv.FV().NumPy()[:] = 0.0
        self.rv.FV().NumPy()[self.idx] = b
        self.sv.data = self.inv * self.rv
        return self.sv.FV().NumPy()[self.idx].copy()

    def solve(self, B, refine=1):
        B = np.asarray(B)
        single = B.ndim == 1
        if single:
            B = B[:, None]
        X = np.zeros(B.shape, dtype=complex)
        for j in range(B.shape[1]):
            b = B[:, j]
            xk = self._apply(b)
            for _ in range(refine):
                self.rv.FV().NumPy()[:] = 0.0
                self.rv.FV().NumPy()[self.idx] = xk
                self.sv.data = self.amat * self.rv
                res = b - self.sv.FV().NumPy()[self.idx]
                nb = np.linalg.norm(b)
                if nb > 0:
                    self.maxres = max(self.maxres, np.linalg.norm(res) / nb)
                xk = xk + self._apply(res)
            X[:, j] = xk
        return X[:, 0] if single else X


def _arc_mid(pr, side):
    w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = pr.p
    if side == "L":
        C, S, E, rad = pr.TL, pr.A1, pr.A2, pL
    else:
        C, S, E, rad = pr.TR, pr.B2, pr.B1, pR
    e1 = np.array(S) - np.array(C)
    e2 = np.array(E) - np.array(C)
    e1 /= np.linalg.norm(e1)
    e2 /= np.linalg.norm(e2)
    bis = e1 + e2
    bis /= np.linalg.norm(bis)
    th = math.acos(max(-1, min(1, float(e1 @ e2))))
    O = np.array(C) + rad / math.sin(th / 2) * bis
    m = np.array(C) - O
    m /= np.linalg.norm(m)
    q = O + rad * m
    return (float(q[0]), float(q[1]))


def _profile_edges(pr):
    from netgen.occ import Pnt, Segment, ArcOfCircle
    zf = 0.5 * (pr.cL + pr.cR)
    edges, corners = [], []
    cur = [(-P / 2, zf)]

    def P2(q):
        return Pnt(q[0], q[1], 0)

    def seg_to(q):
        if math.hypot(q[0] - cur[0][0], q[1] - cur[0][1]) > 1e-9:
            edges.append(Segment(P2(cur[0]), P2(q)))
            cur[0] = q

    def arc_to(mid, q):
        edges.append(ArcOfCircle(P2(cur[0]), P2(mid), P2(q)))
        cur[0] = q

    pv = pr.p
    seg_to(pr.FL)
    corners.append(pr.FL)
    seg_to(pr.A1)
    if pr.dL > 0:
        if pv[8] > 0:
            arc_to(_arc_mid(pr, "L"), pr.A2)
        else:
            corners += [pr.A1, pr.A2]
            seg_to(pr.A2)
    else:
        corners.append(pr.TL)
    seg_to(pr.B2)
    if pr.dR > 0:
        if pv[9] > 0:
            arc_to(_arc_mid(pr, "R"), pr.B1)
        else:
            corners += [pr.B2, pr.B1]
            seg_to(pr.B1)
    else:
        corners.append(pr.TR)
    seg_to(pr.FR)
    corners.append(pr.FR)
    seg_to((P / 2, zf))
    return edges, corners


def build_geometry(p, zb=-300.0, maxh_si=40.0, maxh_ox=80.0, split0=True):
    """Cell [-P/2, P/2] x [zb, H] split (general fuse) by the profile and (if split0) the line
    z = 0."""
    from netgen.occ import Pnt, Segment, Wire, Face, Glue, IdentificationType, gp_Trsf, gp_Vec
    pr = G.Profile(p)
    H = pr.h + pr.t
    edges_prof, corners = _profile_edges(pr)
    xl, xr = -P / 2, P / 2

    def P2(x, z):
        return Pnt(x, z, 0)

    box = Face(Wire([Segment(P2(xl, zb), P2(xr, zb)), Segment(P2(xr, zb), P2(xr, H)),
                     Segment(P2(xr, H), P2(xl, H)), Segment(P2(xl, H), P2(xl, zb))]))
    shape = Glue([box] + edges_prof + ([Segment(P2(xl, 0.0), P2(xr, 0.0))] if split0 else []))
    for f in shape.faces:
        bb = f.bounding_box
        zmax = bb[1][1]
        if split0 and zmax <= 1e-5:
            f.name = "si_sub"
            f.maxh = maxh_si
        elif zmax >= H - 1e-5:
            f.name = "oxide"
            f.maxh = maxh_ox
        else:
            f.name = "si_top"
            f.maxh = maxh_si
    for e in shape.edges:
        if abs(e.start.y - H) < 1e-6 and abs(e.end.y - H) < 1e-6:
            e.name = "top"
        if abs(e.start.y - zb) < 1e-6 and abs(e.end.y - zb) < 1e-6:
            e.name = "bot"
    trf = gp_Trsf.Translation(gp_Vec(P, 0, 0))
    left = [e for e in shape.edges if abs(e.start.x - xl) < 1e-6 and abs(e.end.x - xl) < 1e-6]
    right = [e for e in shape.edges if abs(e.start.x - xr) < 1e-6 and abs(e.end.x - xr) < 1e-6]
    for el in left:
        zs = sorted([el.start.y, el.end.y])
        for er in right:
            zr = sorted([er.start.y, er.end.y])
            if abs(zs[0] - zr[0]) < 1e-6 and abs(zs[1] - zr[1]) < 1e-6:
                el.Identify(er, "periodic", IdentificationType.PERIODIC, trf)
    for vtx in shape.vertices:
        for (cx, cz) in corners:
            if abs(vtx.p.x - cx) < 1e-7 and abs(vtx.p.y - cz) < 1e-7:
                vtx.hpref = 1
    return shape, dict(H=H, zb=zb, corners=corners, profile=pr)


def _mesh_ok(mesh, info, split0):
    """Area check: the mesh must fill the cell and the silicon region exactly (catches netgen
    runs that silently leave faces unmeshed)."""
    import ngsolve as ng
    pr = info["profile"]
    H, zb = info["H"], info["zb"]
    tot = ng.Integrate(1, mesh)
    if abs(tot - P * (H - zb)) > 1e-6 * P * (H - zb):
        return False
    si = sum(ng.Integrate(1, mesh, definedon=mesh.Materials(m)) for m in set(mesh.GetMaterials()) if m.startswith("si"))
    # exact silicon area: P*(-zb) below 0, plus the area under the profile above 0
    _, th, xl, xr, rep = G.layers(pr.p, dz_max=0.5, dx_max=1.0)
    above = float(np.sum(th * (xr - xl) * rep))
    zmin = min(pr.cL, pr.cR)
    exact = P * (-zb) + P * zmin + above
    return abs(si - exact) <= 1e-3 * exact


def make_mesh(p, order=8, maxh_si=40.0, maxh_ox=80.0, hp_levels=4, hp_factor=0.15, zb=-300.0, split0=True):
    import ngsolve as ng
    from netgen.occ import OCCGeometry
    last = None
    for attempt in range(8):
        f = 1.0 - 0.037 * attempt
        try:
            shape, info = build_geometry(p, zb=zb, maxh_si=maxh_si * f, maxh_ox=maxh_ox * f, split0=split0)
            occ = OCCGeometry(shape, dim=2)
            mesh = ng.Mesh(occ.GenerateMesh(maxh=maxh_ox * f))
            if not _mesh_ok(mesh, info, split0):
                last = "area check failed"
                continue
            break
        except Exception as e:          # netgen occasionally fails on near-degenerate inputs
            last = e
    else:
        raise RuntimeError("meshing failed for %r: %r" % (list(p), last))
    if hp_levels > 0:
        mesh.RefineHP(levels=hp_levels, factor=hp_factor)
    mesh.Curve(order)
    return mesh, info


@_threaded
def solve_rg(p, channels, m_inc=7, z_top=None, order=8, maxh_si=40.0, maxh_ox=80.0, hp_levels=4,
             hp_factor=0.15, zb=-300.0, n_top=None, n_bot=40, verbose=False):
    """Reflection matrix Rg[n, m] (|n|, |m| <= m_inc) of the grating (silicon below the profile,
    oxide half-space above, no cover) referenced at the plane z = z_top (default apex + 100 nm):
    for a downward wave exp(i kap m x - i kz_m (z - z_top)) in oxide the upward reflected
    amplitudes Rg[n, m] of exp(i kap n x + i kz_n (z - z_top)).  Total-field formulation with
    exact DtN conditions at z_top (oxide) and zb (silicon)."""
    import ngsolve as ng
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    pr = G.Profile(p)
    if z_top is None:
        # DtN plane 100 nm above the apex: the evanescent orders the corners excite have decayed
        # enough for m_inc + 25 orders; a plane closer to the apex needs orders the boundary mesh
        # cannot resolve
        z_top = pr.apex + 100.0
    q = np.array(p, dtype=float).copy()
    # build the mesh with the cell top at z_top: use t' so that h + t' = z_top
    q[2] = z_top - q[1]
    mesh, info = make_mesh(q, order=order, maxh_si=maxh_si, maxh_ox=maxh_ox, hp_levels=hp_levels,
                           hp_factor=hp_factor, zb=zb, split0=False)
    if n_top is None:
        n_top = m_inc + 25
    fes = ng.Periodic(ng.H1(mesh, order=order, complex=True))
    u, v = fes.TnT()
    x, z = ng.x, ng.y
    free = np.array(list(fes.FreeDofs()), dtype=bool)
    idx = np.nonzero(free)[0]
    kap = 2 * math.pi / P
    nt = np.arange(-n_top, n_top + 1)
    nb = np.arange(-n_bot, n_bot + 1)

    def boundary_modes(name, ns):
        cols = []
        for n in ns:
            lf = ng.LinearForm(fes)
            lf += ng.exp(-1j * kap * n * x) * v * ng.ds(name)
            lf.Assemble()
            cols.append(lf.vec.FV().NumPy()[idx].copy())
        return np.array(cols).T

    Bt = boundary_modes("top", nt)
    Bb = boundary_modes("bot", nb)
    Rgt = np.hstack([Bt, Bb])
    Lft = np.hstack([Bt[:, ::-1], Bb[:, ::-1]])
    ms = np.arange(-m_inc, m_inc + 1)
    out = []
    for (lam, pol) in channels:
        k0 = 2 * math.pi / lam
        esi, eox = complex(eps_si(lam)), complex(eps_sio2(lam))
        kox = k0 * np.sqrt(eox)
        ksi = complex(k0 * np.sqrt(esi))
        if ksi.imag < 0:
            ksi = -ksi
        eps = mesh.MaterialCF({"si_sub": esi, "si_top": esi, "oxide": eox}, default=eox)
        bsi = 1.0 if pol == "TE" else 1.0 / esi
        box = 1.0 if pol == "TE" else 1.0 / eox
        a = ng.BilinearForm(fes, symmetric=True)
        if pol == "TE":
            a += (ng.grad(u) * ng.grad(v) - k0 ** 2 * eps * u * v) * ng.dx
        else:
            a += ((1.0 / eps) * ng.grad(u) * ng.grad(v) - k0 ** 2 * u * v) * ng.dx
        a += (-1j * box * kox * u * v) * ng.ds("top")
        a += (-1j * bsi * ksi * u * v) * ng.ds("bot")
        a.Assemble()
        kzt = np.sqrt(kox ** 2 - (kap * nt) ** 2 + 0j)
        kzt = np.where(kzt.imag < 0, -kzt, kzt)
        kzb = np.sqrt(ksi ** 2 - (kap * nb) ** 2 + 0j)
        kzb = np.where(kzb.imag < 0, -kzb, kzb)
        D = np.concatenate([-(1j / P) * box * (kzt - kox), -(1j / P) * bsi * (kzb - ksi)])
        keep = np.abs(D) > 0
        # right-hand sides: boundary term int beta g v, g = -2 i kz_m exp(i kap m x) on top
        # (total field u = u_inc + u_out, u_out outgoing).  int exp(i kap m x) v ds = l_{-m}(v).
        F = np.zeros((len(idx), len(ms)), dtype=complex)
        kzm = np.sqrt(kox ** 2 - (kap * ms) ** 2 + 0j)
        kzm = np.where(kzm.imag < 0, -kzm, kzm)
        for j, m in enumerate(ms):
            col = np.where(nt == -m)[0][0]
            F[:, j] = box * (-2j * kzm[j]) * Bt[:, col]
        lu = _Solver(a.mat, fes, idx)
        KF = lu.solve(F)
        KL = lu.solve(Lft[:, keep])
        Cm = np.diag(1.0 / D[keep]) + Rgt[:, keep].T @ KL
        Y = np.linalg.solve(Cm, Rgt[:, keep].T @ KF)
        U = KF - KL @ Y
        # reflected amplitudes at z_top: a_n = l_n(u)/P - delta_nm (incident part)
        Rg = np.zeros((len(ms), len(ms)), dtype=complex)
        for i, n in enumerate(ms):
            col = np.where(nt == n)[0][0]
            Rg[i, :] = (Bt[:, col] @ U) / P
        Rg -= np.eye(len(ms))
        out.append(dict(lam=lam, pol=pol, Rg=Rg, kz=kzm, ms=ms))
        if verbose:
            print(lam, pol, "ndof", len(idx), flush=True)
    return dict(z_top=z_top, channels=out, ndof=len(idx))


def cover(Rg_entry, d_ox, lam, pol):
    """+-1 amplitudes in air (referenced at the oxide plane) for a grating with reflection matrix
    Rg (seen from oxide, referenced at a plane at distance d_ox below the oxide plane), for unit
    normal incidence from air."""
    Rg, ms = Rg_entry["Rg"], Rg_entry["ms"]
    k0 = 2 * math.pi / lam
    eox = complex(eps_sio2(lam))
    kap = 2 * math.pi / P
    kzo = np.sqrt(k0 ** 2 * eox - (kap * ms) ** 2 + 0j)
    kzo = np.where(kzo.imag < 0, -kzo, kzo)
    kza = np.sqrt(k0 ** 2 - (kap * ms) ** 2 + 0j)
    kza = np.where(kza.imag < 0, -kza, kza)
    bo = 1.0 if pol == "TE" else 1.0 / eox
    # oxide -> air at the plane: reflection rho (back into oxide), transmission tau_up;
    # air -> oxide transmission tau_dn
    rho = (bo * kzo - kza) / (bo * kzo + kza)
    tau_up = 1 + rho
    tau_dn = 2 * kza / (kza + bo * kzo)
    Phi = np.diag(np.exp(1j * kzo * d_ox))
    n0 = int(np.where(ms == 0)[0][0])
    e0 = np.zeros(len(ms), dtype=complex)
    e0[n0] = tau_dn[n0]
    b = np.linalg.solve(np.eye(len(ms)) - Phi @ np.diag(rho) @ Phi @ Rg, Phi @ e0)
    a_top = Phi @ (Rg @ b)
    A = tau_up * a_top
    ip = int(np.where(ms == 1)[0][0])
    im = int(np.where(ms == -1)[0][0])
    return A[ip], A[im]


# --------------------------------------------------------------------------------------------
# What a channel measures (MODEL.md), at one of two resolutions.
# --------------------------------------------------------------------------------------------
LEVELS = {8: dict(order=8, maxh_si=40.0, maxh_ox=80.0, hp_levels=5, n_bot=60, m_inc=9),
          10: dict(order=10, maxh_si=25.0, maxh_ox=50.0, hp_levels=6, n_bot=80, m_inc=9)}


def amplitudes(marks, channels, order=8):
    """A+, A-, e+, e- (arrays n_marks x n_channels) for marks given as 12 parameters each."""
    marks = np.atleast_2d(np.asarray(marks, dtype=float))
    nm, nc = len(marks), len(channels)
    Ap = np.zeros((nm, nc), dtype=complex)
    Am = np.zeros((nm, nc), dtype=complex)
    for i, p in enumerate(marks):
        rg = solve_rg(p, channels, **LEVELS[order])
        d_ox = p[1] + p[2] - rg["z_top"]
        for j, (lam, pol) in enumerate(channels):
            Ap[i, j], Am[i, j] = cover(rg["channels"][j], d_ox, lam, pol)
    k0 = 2 * math.pi / np.array([c[0] for c in channels])
    kz1 = np.sqrt(k0 ** 2 - (2 * math.pi / P) ** 2)
    ep = np.abs(Ap) ** 2 * (kz1 / k0)[None, :]
    em = np.abs(Am) ** 2 * (kz1 / k0)[None, :]
    return dict(Ap=Ap, Am=Am, ep=ep, em=em)
