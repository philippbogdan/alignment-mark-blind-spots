"""Periodic Helmholtz Green's functions (period P, Bloch wavenumber 0) for complex k.

G(d) = (i/4) sum_m H0^(1)(k |d - m P e_x|),   (Delta + k^2) G = -sum_m delta(d - m P e_x).

Ewald splitting with a frequency-adapted parameter E (Linton 1998, Capolino et al. 2005):

  G = G_spat + G_spec
  G_spat(d) = (1/4pi) sum_{q>=0} c^q/q! E_{q+1}(r^2 E^2),  c = k^2/(4E^2)   (only the image m=0
            after reducing dx to [-P/2, P/2]; E P >= 13.4 makes all other images < 1e-20)
  G_spec(d) = (1/4P) sum_n exp(i kap_n dx) F_n(dz)/gamma_n,
  F_n(z) = exp(gamma z) erfc(a + zE) + exp(-gamma z) erfc(a - zE),  a = gamma/(2E),
  gamma_n = -i kz_n,  kz_n = sqrt(k^2 - kap_n^2) with Im kz_n >= 0.

E is chosen so that the cancellation between both parts (factor exp(Re k^2 / 4E^2)) stays below
`loss` (default 30, i.e. ~1.5 digits).  The z-dependence F_n/gamma_n (and F_n'/gamma_n) is
tabulated once per medium as a Chebyshev series in dz on [-Zr, Zr]; evaluation of the spectral
part for many pairs then costs two small matrix products.

All quantities are returned as functions of d = x - y (target minus source).
"""
import math
import numpy as np
from scipy import special as sp
from scipy.fft import dct

EULER = 0.57721566490153286061


def kz_of(k, kap):
    """sqrt(k^2 - kap^2) on the branch Im >= 0 (outgoing / decaying)."""
    v = np.asarray(k * k - kap * kap, dtype=complex)
    kz = np.sqrt(v)
    kz = np.where(kz.imag < 0, -kz, kz)
    # lossless and propagating: positive real
    kz = np.where((kz.imag == 0) & (kz.real < 0), -kz, kz)
    return kz


class EwaldGreen:
    def __init__(self, k, P, Zr, loss=30.0, tol=1e-17, Emin_xP=13.4, cheb_max=400, cheb_tol=3e-15):
        self.k = complex(k)
        self.P = float(P)
        self.Zr = float(Zr)
        ak = abs(self.k)
        E = max(ak / (2.0 * math.sqrt(math.log(loss))), Emin_xP / P)
        self.E = E
        c = self.k ** 2 / (4 * E * E)
        self.c = c
        # spatial q-series coefficients t_q = c^q / q!
        big = math.exp(max(abs(c), 0.0))
        t = [1.0 + 0j]
        q = 0
        while True:
            q += 1
            t.append(t[-1] * c / q)
            if abs(t[-1]) < tol * 1e-2 * big and q > 4:
                break
        self.tq = np.array(t)
        self.Q = len(t) - 1
        # spectral orders: need kap_N/(2E) >= sqrt(log(1/tol) + Re c) (+ margin)
        amin = math.sqrt(math.log(1.0 / tol) + max(c.real, 0.0)) + 0.5
        N = int(math.ceil(2 * E * amin * P / (2 * math.pi))) + 1
        self.N = N
        n = np.arange(N + 1)
        self.kap = 2 * math.pi * n / P
        kz = kz_of(self.k, self.kap)
        self.kz = kz
        gam = -1j * kz
        self.gam = gam
        self.eps_n = np.where(n == 0, 1.0, 2.0) / (4 * P)
        a = gam / (2 * E)
        self.ea2 = np.exp(-a * a)
        # Chebyshev tables of Phi0 = F/gamma (even) and Phi1 = F'/gamma (odd) on dz = Zr s
        M = 2 * cheb_max
        j = np.arange(M)
        s = np.cos(math.pi * (j + 0.5) / M)
        z = Zr * s
        Ez = (z * E)[:, None]
        wp = a[None, :] + Ez
        wm = a[None, :] - Ez
        pref = np.exp(-(a * a)[None, :] - (z * z * E * E)[:, None])
        xp, xm = sp.erfcx(wp), sp.erfcx(wm)
        Phi0 = pref * (xp + xm) / gam[None, :]
        Phi1 = pref * (xp - xm)
        c0 = dct(Phi0, type=2, axis=0) / M
        c1 = dct(Phi1, type=2, axis=0) / M
        c0[0] /= 2
        c1[0] /= 2
        c0[1::2] = 0.0   # parity (exact)
        c1[0::2] = 0.0
        colscale = np.maximum(np.abs(c0).max(axis=0), np.abs(c1).max(axis=0))[None, :]
        a0 = np.abs(c0) / colscale
        a1 = np.abs(c1) / colscale
        mag = np.maximum(a0.max(axis=1), a1.max(axis=1))
        scale = 1.0
        self._mag = mag
        nz = np.nonzero(mag > cheb_tol)[0]
        L = int(nz[-1]) + 1
        if L > M // 2 - 4 and cheb_tol > 0:
            raise RuntimeError("Chebyshev table for spectral part not converged (L=%d)" % L)
        self.L = L
        self.C0e = np.ascontiguousarray(c0[0:L:2])
        self.C1o = np.ascontiguousarray(c1[1:L:2])
        self.cheb_tail = float(mag[L - 1] / scale) if L > 0 else 0.0
        # values at d = 0 (for diagonal limits)
        self.spec0 = complex(np.sum(self.eps_n * self._phi0_at0()))
        self.spec_xx0 = complex(np.sum(self.eps_n * (-self.kap ** 2) * self._phi0_at0()))
        self.spec_zz0 = complex(np.sum(self.eps_n * (gam ** 2 * self._phi0_at0()
                                                       - 4 * E / math.sqrt(math.pi) * self.ea2)))
        qq = np.arange(1, self.Q + 1)
        # M(0) = lim_{r->0} [G_spat + (1/2pi) log r]
        self.spat_reg0 = (-EULER - 2 * math.log(E) + np.sum(self.tq[1:] / qq)) / (4 * math.pi)
        q2 = np.arange(2, self.Q + 1)
        self.sum_q2 = np.sum(self.tq[2:] / (q2 - 1.0))

    def _phi0_at0(self):
        a = self.gam / (2 * self.E)
        return 2 * self.ea2 * sp.erfcx(a) / self.gam

    # ------------------------------------------------------------------ spectral part
    def spectral(self, dx, dz, hess=True):
        """Returns G, Gx, Gz, Gxx, Gxz, Gzz of the spectral part (1-D arrays)."""
        s = np.clip(dz / self.Zr, -1.0, 1.0)
        L = self.L
        T = np.empty((s.size, L))
        T[:, 0] = 1.0
        if L > 1:
            T[:, 1] = s
        for l in range(2, L):
            T[:, l] = 2 * s * T[:, l - 1] - T[:, l - 2]
        P0 = T[:, 0:L:2] @ self.C0e          # (np, N+1)
        P1 = T[:, 1:L:2] @ self.C1o
        arg = np.outer(dx, self.kap)
        cs, sn = np.cos(arg), np.sin(arg)
        en = self.eps_n
        G = (cs * P0) @ en
        Gz = (cs * P1) @ en
        Gx = -(sn * P0) @ (en * self.kap)
        if not hess:
            return G, Gx, Gz
        Gxx = -(cs * P0) @ (en * self.kap ** 2)
        Gxz = -(sn * P1) @ (en * self.kap)
        Gzz = (cs * P0) @ (en * self.gam ** 2) \
            - (4 * self.E / math.sqrt(math.pi)) * np.exp(-(dz * self.E) ** 2) * (cs @ (en * self.ea2))
        return G, Gx, Gz, Gxx, Gxz, Gzz

    # ------------------------------------------------------------------ spatial part
    def spatial(self, dxr, dz, r2):
        """Spatial (m=0) part.  Returns dict with
        g0  = (1/4pi) E_1(x)                      (q = 0 term, k-independent)
        gq  = (1/4pi) sum_{q>=1} t_q E_{q+1}(x)
        a1  : grad of full spatial part = a1 * d
        h1, h2 : Hess of q>=1 part = h1 I + h2 d d^T
        emx = exp(-x)
        """
        E2 = self.E ** 2
        x = r2 * E2
        cut = x < 80.0
        xs = np.where(cut, x, 80.0)
        emx = np.exp(-xs)
        E1 = sp.exp1(xs)
        E0 = emx / xs
        # recurrence E_{n+1} = (e^{-x} - x E_n)/n
        Em1 = E0          # E_{q-1} for q = 1
        Eq = E1           # E_q    for q = 1
        S0 = np.zeros(x.shape, complex)   # sum t_q E_{q+1}
        S1 = np.zeros(x.shape, complex)   # sum t_q E_q
        S2 = np.zeros(x.shape, complex)   # sum t_q E_{q-1}
        for q in range(1, self.Q + 1):
            Eq1 = (emx - xs * Eq) / q
            tq = self.tq[q]
            S0 += tq * Eq1
            S1 += tq * Eq
            S2 += tq * Em1
            Em1, Eq = Eq, Eq1
        inv4pi = 1.0 / (4 * math.pi)
        g0 = np.where(cut, inv4pi * E1, 0.0)
        gq = np.where(cut, inv4pi * S0, 0.0)
        # gradient: q=0: -(1/2pi) e^{-x}/r^2 ; q>=1: -(E^2/2pi) S1
        a1 = np.where(cut, -emx / (2 * math.pi * r2) - (E2 / (2 * math.pi)) * S1, 0.0)
        h1 = np.where(cut, -(E2 / (2 * math.pi)) * S1, 0.0)
        h2 = np.where(cut, (E2 * E2 / math.pi) * S2, 0.0)
        return dict(g0=g0, gq=gq, a1=a1, h1=h1, h2=h2, emx=np.where(cut, emx, 0.0), x=x)

    def reduce_dx(self, dx):
        return dx - self.P * np.round(dx / self.P)

    def full(self, dx, dz):
        """G, Gx, Gz (full periodic) -- for testing."""
        dxr = self.reduce_dx(dx)
        r2 = dxr ** 2 + dz ** 2
        sp_ = self.spatial(dxr, dz, r2)
        G, Gx, Gz = self.spectral(dx, dz, hess=False)
        return (G + sp_["g0"] + sp_["gq"], Gx + sp_["a1"] * dxr, Gz + sp_["a1"] * dz)


def q0_hess_diff(g1, g2, dxr, dz, r2):
    """Hessian of [(1/4pi)E_1(r^2 E1^2) - (1/4pi)E_1(r^2 E2^2)], stable for r -> 0.
    Returns (Hxx, Hxz, Hzz)."""
    E1s, E2s = g1.E ** 2, g2.E ** 2
    x2 = r2 * E2s
    e2 = np.exp(-np.minimum(x2, 700.0))
    D1 = e2 * np.expm1(-(E1s - E2s) * r2) / r2           # (e^{-x1} - e^{-x2}) / r^2
    D2 = E1s * np.exp(-np.minimum(r2 * E1s, 700.0)) - E2s * e2
    ux, uz = dxr * dxr / r2, dz * dz / r2
    uxz = dxr * dz / r2
    f = -1.0 / (2 * math.pi)
    Hxx = f * (D1 * (1 - 2 * ux) - 2 * D2 * ux)
    Hzz = f * (D1 * (1 - 2 * uz) - 2 * D2 * uz)
    Hxz = f * (D1 * (-2 * uxz) - 2 * D2 * uxz)
    return Hxx, Hxz, Hzz


class LayerReflection:
    """Smooth reflected part of the oxide Green's function from the oxide/air plane z = H:
    G_ref = sum_n (i/(2 P kz_n)) R_n exp(i kap_n (x - x')) exp(i kz_n (2H - z - z')),
    R_n = (b1 kz1 - b0 kz0)/(b1 kz1 + b0 kz0) (oxide -> air, beta-weighted)."""

    def __init__(self, k_ox, k_air, b_ox, b_air, H, P, zmax, tol=1e-18):
        self.H, self.P = H, P
        Dmin = 2 * (H - zmax)
        assert Dmin > 0
        n = 0
        while True:
            n += 1
            kap = 2 * math.pi * n / P
            kz = kz_of(k_ox, kap)
            if abs(np.exp(1j * kz * Dmin)) < tol and kap > abs(k_ox):
                break
        N = n
        nn = np.arange(-N, N + 1)
        self.kap = 2 * math.pi * nn / P
        kz1 = kz_of(k_ox, self.kap)
        kz0 = kz_of(k_air, self.kap)
        R = (b_ox * kz1 - b_air * kz0) / (b_ox * kz1 + b_air * kz0)
        self.kz = kz1
        self.coef = 1j / (2 * P * kz1) * R

    def mats(self, xt, zt, xs, zs, hess=True):
        """Return G, G_X, G_Z, (G_XX, G_XZ, G_ZZ) with X = x - x', Z = z + z' as (nt, ns) arrays."""
        Ut = np.exp(1j * np.outer(xt, self.kap) + 1j * np.outer(self.H - zt, self.kz))
        Us = np.exp(-1j * np.outer(xs, self.kap) + 1j * np.outer(self.H - zs, self.kz))
        c = self.coef
        ik = 1j * self.kap
        mz = -1j * self.kz
        out = [Ut @ (c[:, None] * Us.T), Ut @ ((c * ik)[:, None] * Us.T), Ut @ ((c * mz)[:, None] * Us.T)]
        if hess:
            out += [Ut @ ((c * ik * ik)[:, None] * Us.T), Ut @ ((c * ik * mz)[:, None] * Us.T),
                    Ut @ ((c * mz * mz)[:, None] * Us.T)]
        return out
