"""Quick unit checks of the BIE building blocks (about 30 s on a laptop).

  python bie/selftest.py

1. Ewald periodic Green's function (value, gradient, Hessian) vs the brute-force spectral sum
   (60001 orders) for |dz| >= 0.6 nm, Si and SiO2, lossy and lossless, 400..1100 nm.
2. Self-panel split K = A log r + B: B(x, y) -> diag_B as y -> x along a curved arc.
3. Small-|d| kernel path vs full Ewald path.
4. Log-product weights vs mpmath.
5. Flat Si/SiO2 interface: reflection vs the analytic thin-film formula; +-1 orders vanish.
"""
import os
import sys
import math
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bie.greens import EwaldGreen, kz_of                 # noqa: E402
from bie.kernels import MullerKernel                     # noqa: E402
from bie import materials, solver                        # noqa: E402
from bie.geometry import Line                            # noqa: E402
from bie.quad import Tables, log_moments                 # noqa: E402

P = 2000.0
worst = {}


def check(name, val, tol):
    worst[name] = val
    print("%-55s %.1e  (tol %.0e) %s" % (name, val, tol, "ok" if val < tol else "FAIL"))


def direct(k, dx, dz, N=30000):
    out = [np.zeros(dx.size, complex) for _ in range(6)]
    for a in range(-N, N + 1, 5000):
        n = np.arange(a, min(a + 5000, N + 1))
        kap = 2 * np.pi * n / P
        kz = kz_of(k, kap)
        sg = np.sign(dz)[:, None]
        ph = (1j / (2 * P)) * np.exp(1j * np.outer(dx, kap) + 1j * np.outer(np.abs(dz), kz)) / kz
        for o, f in zip(out, [1, 1j * kap, 1j * kz * sg, (1j * kap) ** 2, (1j * kap) * (1j * kz * sg), (1j * kz) ** 2]):
            o += (ph * f).sum(1)
    return out


def ewald_all(g, dx, dz):
    dxr = g.reduce_dx(dx)
    r2 = dxr ** 2 + dz ** 2
    s = g.spatial(dxr, dz, r2)
    G, Gx, Gz, Gxx, Gxz, Gzz = g.spectral(dx, dz)
    emx, E2 = s["emx"], g.E ** 2
    f = -1 / (2 * np.pi) * emx
    q0xx = f * (1 / r2 - 2 * dxr * dxr / r2 ** 2 - 2 * E2 * dxr * dxr / r2)
    q0zz = f * (1 / r2 - 2 * dz * dz / r2 ** 2 - 2 * E2 * dz * dz / r2)
    q0xz = f * (-2 * dxr * dz / r2 ** 2 - 2 * E2 * dxr * dz / r2)
    return [G + s["g0"] + s["gq"], Gx + s["a1"] * dxr, Gz + s["a1"] * dz,
            Gxx + q0xx + s["h1"] + s["h2"] * dxr * dxr, Gxz + q0xz + s["h2"] * dxr * dz,
            Gzz + q0zz + s["h1"] + s["h2"] * dz * dz]


def test_green():
    rng = np.random.default_rng(0)
    dx = np.concatenate([rng.uniform(-2500, 2500, 40), [0.3, -2.0, 10.0, 999.0, 1999.7]])
    dz = np.concatenate([rng.uniform(-290, 290, 40), [0.6, 1.0, -2.0, 3.0, 1.0]])
    dz = np.where(np.abs(dz) < 0.6, 0.6, dz)
    err = 0.0
    for lam in (400.0, 633.0, 1100.0):
        for k in (2 * np.pi * materials.n_si(lam) / lam, 2 * np.pi * materials.n_si(lam).real / lam,
                  2 * np.pi * materials.n_sio2(lam) / lam):
            g = EwaldGreen(k, P, 300.0)
            A, B = ewald_all(g, dx, dz), direct(k, dx, dz)
            for a, b in zip(A, B):
                err = max(err, np.max(np.abs(a - b) / (np.abs(b) + np.abs(b).max() * 1e-3)))
    check("Ewald G, grad G, Hess G vs brute-force sum (rel)", err, 1e-12)


def make_kernel(lam, pol):
    k0 = 2 * np.pi / lam
    e1, e2 = materials.eps_sio2(lam), materials.eps_si(lam)
    b1, b2 = (1, 1) if pol == "TE" else (1 / e1, 1 / e2)
    return MullerKernel(k0 * np.sqrt(e1), k0 * np.sqrt(e2), k0, b1, b2, 1.0, 960.0, P, 0.0, 260.0)


def test_split():
    err = 0.0
    for pol in ("TE", "TM"):
        K = make_kernel(633.0, pol)
        R, sg, c, phi0 = 50.0, -1.0, np.array([10.0, 100.0]), 0.7

        def pt(phi):
            x, z = c[0] + R * np.cos(phi), c[1] + R * np.sin(phi)
            tx, tz = -sg * np.sin(phi), sg * np.cos(phi)
            a = np.atleast_1d
            return dict(rx=a(x), rz=a(z), x=a(x), z=a(z), nx=a(-tz), nz=a(tx))
        T = pt(phi0)
        Bd = K.diag_B(T["nx"], T["nz"], np.array([sg / R]))
        prev = None
        for d in (1e-2, 1e-3):
            S = pt(phi0 + sg * d / R)
            M = K.matrix(T, S, with_ref=False)
            dx, dz = T["rx"] - S["rx"], T["rz"] - S["rz"]
            A = K.log_coeffs(dx, dz, T["nx"], T["nz"], S["nx"], S["nz"])
            r = np.hypot(dx, dz)
            e = [abs((M[i] - A[i] * np.log(r))[0, 0] - Bd[i][0]) for i in range(4)]
            if prev is not None:   # linear convergence in d: error ratio ~ 10
                err = max(err, max(ei / max(pi_, 1e-300) for ei, pi_ in zip(e, prev) if pi_ > 1e-12))
            prev = e
    check("diag limit of B: error ratio for d 1e-2 -> 1e-3 (expect 0.1)", err, 0.15)


def test_small():
    err = 0.0
    for pol in ("TE", "TM"):
        K = make_kernel(400.0, pol)
        rng = np.random.default_rng(1)
        for sc in (K.r_small * 0.9, 1e-6):
            dx, dz = rng.uniform(-sc, sc, 50), rng.uniform(-sc, sc, 50)
            nx, nz, mx, mz = [rng.normal(size=50) for _ in range(4)]
            a = K._periodic_small(dx, dz, nx, nz, mx, mz)
            rs = K.r_small
            K.r_small = 0.0
            b = K._periodic(dx, dz, nx, nz, mx, mz)
            K.r_small = rs
            r = np.hypot(dx, dz)
            err = max(err, max(np.max(np.abs(x - y) * r) for x, y in zip(a, b)))
    check("small-|d| path vs Ewald path (|dK| * r)", err, 1e-13)


def test_quad():
    import mpmath as mp
    mp.mp.dps = 30
    T = Tables()
    err = 0.0
    for a in (T.t[0], T.t[7], 0.3):
        ex = [mp.quad(lambda t: mp.legendre(k, t) * mp.log(abs(t - a)), [-1, a, 1]) for k in range(16)]
        err = max(err, max(abs(float(e) - m) for e, m in zip(ex, log_moments(a, 16))))
    check("Legendre log moments vs mpmath", err, 1e-14)


class Flat:
    def __init__(self, z0=50.0, H=960.0):
        self.segs = [Line((-300.0, z0), (1700.0, z0))]
        self.H, self.zmin, self.zmax = H, z0, z0 + 1e-9


def test_flat():
    err, eo = 0.0, 0.0
    for lam in (450.0, 850.0):
        for pol in ("TE", "TM"):
            pr = Flat()
            r = solver.solve_one(pr, lam, pol, dict(solver.DEFAULTS), want_transmission=True)
            k0 = 2 * np.pi / lam
            e1, e2 = materials.eps_sio2(lam), materials.eps_si(lam)
            k1, k2 = k0 * np.sqrt(e1), k0 * np.sqrt(e2)
            b = (1, 1, 1) if pol == "TE" else (1, 1 / e1, 1 / e2)
            r01 = (b[0] * k0 - b[1] * k1) / (b[0] * k0 + b[1] * k1)
            r12 = (b[1] * k1 - b[2] * k2) / (b[1] * k1 + b[2] * k2)
            ph = np.exp(2j * k1 * (pr.H - pr.zmin))
            rex = (r01 + r12 * ph) / (1 + r01 * r12 * ph)
            err = max(err, abs(r["orders"][0] + r01 - rex))
            eo = max(eo, abs(r["Ap"]), abs(r["Am"]))
    check("flat interface: |r - r_thin_film|", err, 1e-13)
    check("flat interface: |A+-|", eo, 1e-13)


if __name__ == "__main__":
    test_green()
    test_split()
    test_small()
    test_quad()
    test_flat()
    try:
        solver.solve([1000, 160, 800, 0, 0, 0, 0, 0, 0, 0, 0, 0.0], [(727.3477662873167, "TE")])
        print("Wood-anomaly guard                                      FAIL")
    except ValueError:
        print("Wood-anomaly guard (oxide order grazing at 727.3478 nm) ok")
