"""Boundary-integral (Mueller / RCIP / Ewald) solver for the alignment mark.

Entry point:  solve(params, channels, **accuracy) -> dict(Ap, Am, ep, em) with arrays of shape
(n_marks, n_channels).  params: one 12-vector or an (n_marks, 12) array
[w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR]; channels: list of (lambda_nm, 'TE'|'TM').
Amplitudes of the reflected +-1 orders in air (TE: E_y, TM: H_y), unit incident plane wave at
normal incidence, exp(-i omega t), order n ~ exp(+2 pi i n x / P), referenced to z = h + t.

accuracy keywords (defaults in DEFAULTS):
  ppw      panels per wavelength in Si (coarse panel length = lambda / (Re n_Si * ppw))
  hcap     absolute cap on the coarse panel length (nm)
  growth   grading rate of panel sizes away from corners
  nlev     RCIP refinement levels per corner (smallest scale s 2^-nlev)
  order    Gauss-Legendre points per panel (16 only: log weights are tabulated for any n)
  ewald_loss  admitted cancellation factor in the Ewald split
  lossless    set Im eps_Si = 0 (energy test)
"""
import math
import time
import numpy as np
from . import materials
from .geometry import Profile, P as PERIOD
from .kernels import MullerKernel
from .quad import Tables
from .nystrom import assemble
from .rcip import coarse_mesh, cyclic_adjacency, compute_R, _comp_idx
from .greens import kz_of

DEFAULTS = dict(ppw=2.5, hcap=60.0, growth=0.5, nlev=55, order=16, ewald_loss=30.0,
                lossless=False, corner_factor=1.0, aux_levels=9, direct_levels=0)
_TABLES = {}


def _tables(order, aux_levels):
    key = (order, aux_levels)
    if key not in _TABLES:
        _TABLES[key] = Tables(order, n_aux=16, aux_levels=aux_levels)
    return _TABLES[key]


def solve_one(prof, lam, pol, acc, want_transmission=False, return_info=False):
    t0 = time.time()
    k0 = 2 * math.pi / lam
    eps1 = complex(materials.eps_sio2(lam))
    eps2 = complex(materials.eps_si(lam))
    if acc["lossless"]:
        eps2 = complex(eps2.real, 0.0)
    k1, k2 = k0 * np.sqrt(eps1), k0 * np.sqrt(eps2)
    if pol == "TE":
        b0, b1, b2 = 1.0, 1.0, 1.0
    elif pol == "TM":
        b0, b1, b2 = 1.0, 1.0 / eps1, 1.0 / eps2
    else:
        raise ValueError(pol)
    H = prof.H
    nW = np.arange(0, int(PERIOD * abs(k1) / (2 * math.pi)) + 2)
    if np.min(np.abs(kz_of(k1, 2 * math.pi * nW / PERIOD))) * PERIOD < 1e-4:
        raise ValueError("lambda = %r nm is (within ~1e-8 nm) an oxide Rayleigh/Wood wavelength; "
                         "the periodic Green's function is singular there" % lam)
    kern = MullerKernel(k1, k2, k0, b1, b2, b0, H, PERIOD, prof.zmin, prof.zmax,
                        ewald_loss=acc["ewald_loss"])
    tab = _tables(acc["order"], acc["aux_levels"])
    n = tab.n
    hmax = min(acc["hcap"], lam / (np.real(np.sqrt(eps2)) * acc["ppw"]))
    direct = acc["direct_levels"] > 0
    panels, corners = coarse_mesh(prof, hmax, acc["growth"], acc["corner_factor"], acc["direct_levels"])
    npan = len(panels)
    N = n * npan
    skip = set()
    if not direct:
        for c in corners:
            for p in c["panels"]:
                for q in c["panels"]:
                    skip.add((p, q))
    K, ND, w = assemble(kern, tab, panels, "abs", cyclic_adjacency(panels, all_pairs=direct), skip_blocks=skip)
    t_asm = time.time() - t0
    # RCIP
    Rs = []
    for c in (corners if not direct else []):
        idx = _comp_idx(c["panels"], n, N)
        K[np.ix_(idx, idx)] = 0.0
        R = compute_R(kern, tab, c, acc["nlev"])
        Rs.append((idx, R))
        K[:, idx] = K[:, idx] @ R
    t_rcip = time.time() - t0 - t_asm
    # right-hand side: incident field in the oxide (background air/oxide half-spaces)
    tau0 = 2 * b0 * k0 / (b0 * k0 + b1 * k1)
    uinc = tau0 * np.exp(-1j * k1 * (ND["z"] - H))
    dninc = ND["nz"] * (-1j * k1) * uinc
    g = np.concatenate([b1 * uinc / kern.c1, dninc / kern.c2])
    K[np.diag_indices_from(K)] += 1.0
    rt = np.linalg.solve(K, g)
    rho = rt.copy()
    for idx, R in Rs:
        rho[idx] = R @ rt[idx]
    u, v = rho[:N], rho[N:]
    # reflected orders in air (amplitude at z = H)
    out = {}
    nn = np.array([-1, 0, 1]) if not want_transmission else None
    if want_transmission:
        nmax = int(math.floor(PERIOD / lam)) + 1
        nn = np.arange(-nmax, nmax + 1)
    kap = 2 * math.pi * nn / PERIOD
    kz0 = kz_of(k0, kap)
    kz1 = kz_of(k1, kap)
    tn = 2 * b1 * kz1 / (b1 * kz1 + b0 * kz0)
    coef = 1j * tn / (2 * PERIOD * kz1)
    ph = np.exp(-1j * np.outer(kap, ND["x"]) + 1j * np.outer(kz1, H - ND["z"]))
    dn = -1j * np.outer(kap, ND["nx"]) - 1j * np.outer(kz1, ND["nz"])
    A = coef * ((ph * (u[None, :] * dn - v[None, :] / b1)) @ w)
    amp = dict(zip(nn.tolist(), A))
    out["Ap"], out["Am"] = amp[1], amp[-1]
    out["ep"] = abs(out["Ap"]) ** 2 * kz0[nn == 1][0].real / k0
    out["em"] = abs(out["Am"]) ** 2 * kz0[nn == -1][0].real / k0
    if want_transmission:
        r0 = (b0 * k0 - b1 * k1) / (b0 * k0 + b1 * k1)
        rr = A.copy()
        rr[nn == 0] += r0
        Rflux = np.sum((b0 * kz0).real * np.abs(rr) ** 2 * (kz0.imag == 0)) / (b0 * k0)
        nmax2 = int(math.floor(PERIOD * abs(k2) / (2 * math.pi))) + 1
        n2 = np.arange(-nmax2, nmax2 + 1)
        kap2 = 2 * math.pi * n2 / PERIOD
        kz2 = kz_of(k2, kap2)
        zr = prof.zmin
        ph2 = np.exp(-1j * np.outer(kap2, ND["x"]) + 1j * np.outer(kz2, ND["z"] - zr))
        dn2 = -1j * np.outer(kap2, ND["nx"]) + 1j * np.outer(kz2, ND["nz"])
        T = (1j / (2 * PERIOD * kz2)) * ((ph2 * (v[None, :] / b2 - u[None, :] * dn2)) @ w)
        prop = (np.abs(kz2.imag) < 1e-14 * abs(k2))
        Tflux = np.sum((b2 * kz2).real * np.abs(T) ** 2 * prop) / (b0 * k0)
        out.update(R=Rflux, T=Tflux, orders=amp)
    out.update(time=time.time() - t0, t_asm=t_asm, t_rcip=t_rcip, N=N, npan=npan,
               ncorners=len(corners))
    if return_info:
        out["info"] = dict(kern=kern, panels=panels, corners=corners, ND=ND, w=w, u=u, v=v)
    return out


def solve(params, channels, **accuracy):
    acc = dict(DEFAULTS)
    for k in accuracy:
        if k not in DEFAULTS:
            raise TypeError("unknown accuracy keyword %s" % k)
    acc.update(accuracy)
    params = np.atleast_2d(np.asarray(params, float))
    nm, nc = params.shape[0], len(channels)
    Ap = np.empty((nm, nc), complex)
    Am = np.empty((nm, nc), complex)
    ep = np.empty((nm, nc))
    em = np.empty((nm, nc))
    tm = np.empty((nm, nc))
    for i, p in enumerate(params):
        prof = Profile(p)
        for j, (lam, pol) in enumerate(channels):
            r = solve_one(prof, float(lam), pol, acc)
            Ap[i, j], Am[i, j], ep[i, j], em[i, j], tm[i, j] = r["Ap"], r["Am"], r["ep"], r["em"], r["time"]
    return dict(Ap=Ap, Am=Am, ep=ep, em=em, time=tm)


def readings(sol):
    """delta (nm), V, e+ + e- from a solve() result."""
    delta = -(PERIOD / (4 * math.pi)) * np.angle(sol["Ap"] / sol["Am"])
    s = sol["ep"] + sol["em"]
    V = 2 * np.sqrt(sol["ep"] * sol["em"]) / s
    return delta, V, s
