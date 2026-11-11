"""Mueller-type kernels for the transmission problem across Gamma.

Medium 1 = SiO2 above Gamma (Green's function G1 = periodic oxide G + reflected part from the
oxide/air plane), medium 2 = Si below (periodic G2).  Unknowns on Gamma: u (field) and
v = beta d_n u (beta = 1 for TE, 1/eps for TM), n pointing into medium 1.

    c1 u + (beta2 K2 - beta1 K1) u + (S1 - S2) v      = beta1 u_inc
    c2 v - (T1 - T2) u + (K1'/beta1 - K2'/beta2) v    = d_n u_inc
    c1 = (beta1 + beta2)/2,  c2 = (1/beta1 + 1/beta2)/2

with S, K (d/dn_y), K' (d/dn_x), T (d^2/dn_x dn_y) the layer operators.  Rows are divided by
c1, c2, giving (I + K~) rho = g~.  The kernels are evaluated as functions of d = x - y with
relative coordinates so that RCIP meshes down to 1e-17 nm are exact; T1 - T2 is formed
without cancellation (the k-independent q=0 Ewald term is differenced analytically).

For the self-panel product quadrature every kernel is split K = A log|x-y| + B, A, B smooth:
A12 = -(1/2pi)(J0(k1 r) - J0(k2 r)), A11, A22 from J1, A21 from the T-kernel, and the
diagonal limits of B are given in `diag_B`.
"""
import math
import numpy as np
from scipy import special as sp
from .greens import EwaldGreen, LayerReflection, q0_hess_diff, EULER

CHUNK = 20000


class MullerKernel:
    def __init__(self, k_ox, k_si, k_air, beta_ox, beta_si, beta_air, H, P, zmin, zmax,
                 ewald_loss=30.0):
        self.k1, self.k2 = complex(k_ox), complex(k_si)
        self.b1, self.b2 = complex(beta_ox), complex(beta_si)
        self.P = P
        self.H = H
        Zr = (zmax - zmin) * 1.02 + 1.0
        self.g1 = EwaldGreen(self.k1, P, Zr, loss=ewald_loss)
        self.g2 = EwaldGreen(self.k2, P, Zr, loss=ewald_loss)
        self.ref = LayerReflection(self.k1, k_air, beta_ox, beta_air, H, P, zmax)
        self.c1 = (self.b1 + self.b2) / 2
        self.c2 = (1 / self.b1 + 1 / self.b2) / 2
        # small-|d| expansion data: G = -(1/2pi) J0(kr) log r + S0 + d^T Hs d / 2 + O(d^4)
        self.small = []
        for g in (self.g1, self.g2):
            hM = (g.E ** 2 / (2 * math.pi)) * (1 + g.tq[1] * (EULER - 1 + 2 * math.log(g.E)) - g.sum_q2)
            self.small.append(dict(S0=g.spat_reg0 + g.spec0, Hxx=hM + g.spec_xx0, Hzz=hM + g.spec_zz0))
        kmax = max(abs(self.k1), abs(self.k2), self.g1.E, self.g2.E)
        self.r_small = 1e-4 / kmax

    # ------------------------------------------------------------------ main evaluation
    def matrix(self, T, S, with_ref=True):
        """Kernel values M11..M22 (unweighted, rows already divided by c1, c2) for all target
        nodes T (dict rx, rz, x, z, nx, nz) against source nodes S.  Pairs with r = 0 get the
        smooth (reflected) part only."""
        nt, ns = T["rx"].size, S["rx"].size
        dx = (T["rx"][:, None] - S["rx"][None, :]).ravel()
        dz = (T["rz"][:, None] - S["rz"][None, :]).ravel()
        tnx = np.broadcast_to(T["nx"][:, None], (nt, ns)).ravel()
        tnz = np.broadcast_to(T["nz"][:, None], (nt, ns)).ravel()
        snx = np.broadcast_to(S["nx"][None, :], (nt, ns)).ravel()
        snz = np.broadcast_to(S["nz"][None, :], (nt, ns)).ravel()
        out = [np.empty(nt * ns, complex) for _ in range(4)]
        for a in range(0, nt * ns, CHUNK):
            b = min(a + CHUNK, nt * ns)
            r = self._periodic(dx[a:b], dz[a:b], tnx[a:b], tnz[a:b], snx[a:b], snz[a:b])
            for o, v in zip(out, r):
                o[a:b] = v
        out = [o.reshape(nt, ns) for o in out]
        if with_ref:
            G, GX, GZ, GXX, GXZ, GZZ = self.ref.mats(T["x"], T["z"], S["x"], S["z"])
            nx, nz = T["nx"][:, None], T["nz"][:, None]
            mx, mz = S["nx"][None, :], S["nz"][None, :]
            Kr = -mx * GX + mz * GZ              # d/dn_y
            Kpr = nx * GX + nz * GZ              # d/dn_x
            Tr = -nx * mx * GXX + nx * mz * GXZ - nz * mx * GXZ + nz * mz * GZZ
            out[0] += (-self.b1 * Kr) / self.c1
            out[1] += G / self.c1
            out[2] += (-Tr) / self.c2
            out[3] += (Kpr / self.b1) / self.c2
        return out

    def _periodic(self, dx, dz, nx, nz, mx, mz):
        if np.all(np.abs(dx) < self.r_small) and np.all(np.abs(dz) < self.r_small):
            return self._periodic_small(dx, dz, nx, nz, mx, mz)
        g1, g2 = self.g1, self.g2
        dxr = g1.reduce_dx(dx)
        r2 = dxr * dxr + dz * dz
        zero = r2 == 0
        r2 = np.where(zero, 1.0, r2)
        s1 = g1.spatial(dxr, dz, r2)
        s2 = g2.spatial(dxr, dz, r2)
        p1 = g1.spectral(dx, dz)
        p2 = g2.spectral(dx, dz)
        # values and gradients
        G1 = s1["g0"] + s1["gq"] + p1[0]
        G2 = s2["g0"] + s2["gq"] + p2[0]
        G1x, G1z = s1["a1"] * dxr + p1[1], s1["a1"] * dz + p1[2]
        G2x, G2z = s2["a1"] * dxr + p2[1], s2["a1"] * dz + p2[2]
        # Hessian difference
        Hxx, Hxz, Hzz = q0_hess_diff(g1, g2, dxr, dz, r2)
        Hxx = Hxx + (s1["h1"] + s1["h2"] * dxr * dxr) - (s2["h1"] + s2["h2"] * dxr * dxr) + p1[3] - p2[3]
        Hzz = Hzz + (s1["h1"] + s1["h2"] * dz * dz) - (s2["h1"] + s2["h2"] * dz * dz) + p1[5] - p2[5]
        Hxz = Hxz + (s1["h2"] - s2["h2"]) * dxr * dz + p1[4] - p2[4]
        K1 = -(mx * G1x + mz * G1z)
        K2 = -(mx * G2x + mz * G2z)
        K1p = nx * G1x + nz * G1z
        K2p = nx * G2x + nz * G2z
        dT = -(nx * mx * Hxx + (nx * mz + nz * mx) * Hxz + nz * mz * Hzz)
        M11 = (self.b2 * K2 - self.b1 * K1) / self.c1
        M12 = (G1 - G2) / self.c1
        M21 = -dT / self.c2
        M22 = (K1p / self.b1 - K2p / self.b2) / self.c2
        for M in (M11, M12, M21, M22):
            M[zero] = 0.0
        return M11, M12, M21, M22

    def _periodic_small(self, dx, dz, nx, nz, mx, mz):
        """Same as _periodic for |d| < r_small (|k d| < 1e-3): exact log part by series,
        smooth part by its second-order Taylor polynomial at d = 0."""
        r2 = dx * dx + dz * dz
        zero = r2 == 0
        r2 = np.where(zero, 1.0, r2)
        r = np.sqrt(r2)
        lr = np.log(r)
        f = 1.0 / (2 * math.pi)
        J = 8
        res = []
        for k, sm in ((self.k1, self.small[0]), (self.k2, self.small[1])):
            # J0(kr) = sum a_j r^(2j)
            a = [((-1) ** j) * (k * k / 4) ** j / math.factorial(j) ** 2 for j in range(J)]
            J0 = sum(a[j] * r2 ** j for j in range(J))
            dJ0 = sum(2 * j * a[j] * r2 ** j for j in range(1, J)) / r      # d/dr J0(kr)
            # f(r) = -(1/2pi) J0 log r ; f' = -(1/2pi)(J0' log r + J0/r)
            G = -f * J0 * lr + sm["S0"] + 0.5 * (sm["Hxx"] * dx * dx + sm["Hzz"] * dz * dz)
            fp_over_r = -f * (dJ0 * lr / r + J0 / r2)
            Gx = fp_over_r * dx + sm["Hxx"] * dx
            Gz = fp_over_r * dz + sm["Hzz"] * dz
            res.append((G, Gx, Gz, a))
        (G1, G1x, G1z, a1), (G2, G2x, G2z, a2) = res
        # Hessian of -(1/2pi)(J0(k1 r) - J0(k2 r)) log r, series without cancellation
        da = [a1[j] - a2[j] for j in range(J)]
        D0 = sum(da[j] * r2 ** (j - 1) for j in range(1, J))            # dJ0 / r^2
        D1 = sum(2 * j * da[j] * r2 ** (j - 1) for j in range(1, J))    # dJ0' / r
        D2 = sum(2 * j * (2 * j - 1) * da[j] * r2 ** (j - 1) for j in range(1, J))   # dJ0''
        # h(r) = -(1/2pi) dJ0 log r: h'/r = -(1/2pi)(D1 log r + D0), h'' = -(1/2pi)(D2 log r + 2 D1 - D0)
        hpr = -f * (D1 * lr + D0)
        hpp = -f * (D2 * lr + 2 * D1 - D0)
        ux, uz, uxz = dx * dx / r2, dz * dz / r2, dx * dz / r2
        s1, s2 = self.small
        Hxx = hpp * ux + hpr * (1 - ux) + (s1["Hxx"] - s2["Hxx"])
        Hzz = hpp * uz + hpr * (1 - uz) + (s1["Hzz"] - s2["Hzz"])
        Hxz = (hpp - hpr) * uxz
        K1 = -(mx * G1x + mz * G1z)
        K2 = -(mx * G2x + mz * G2z)
        K1p = nx * G1x + nz * G1z
        K2p = nx * G2x + nz * G2z
        dT = -(nx * mx * Hxx + (nx * mz + nz * mx) * Hxz + nz * mz * Hzz)
        M11 = (self.b2 * K2 - self.b1 * K1) / self.c1
        M12 = (G1 - G2) / self.c1
        M21 = -dT / self.c2
        M22 = (K1p / self.b1 - K2p / self.b2) / self.c2
        for M in (M11, M12, M21, M22):
            M[zero] = 0.0
        return M11, M12, M21, M22

    # ------------------------------------------------------------------ log-split pieces
    def log_coeffs(self, dx, dz, nx, nz, mx, mz):
        """A-coefficients of log|x-y| for the four kernels (same shapes as inputs, r > 0 or r = 0)."""
        r = np.hypot(dx, dz)
        rs = np.where(r == 0, 1.0, r)
        ux, uz = dx / rs, dz / rs
        nd = np.where(r == 0, 0.0, nx * ux + nz * uz)     # n_x . d^
        md = np.where(r == 0, 0.0, mx * ux + mz * uz)     # n_y . d^
        nm = nx * mx + nz * mz
        out = []
        J = {}
        for key, k in (("1", self.k1), ("2", self.k2)):
            z = k * r
            J0 = sp.jv(0, z)
            J1 = sp.jv(1, z)
            zs = np.where(z == 0, 1.0, z)
            J1oz = np.where(z == 0, 0.5, J1 / zs)          # J1(z)/z
            J[key] = (J0, J1oz)
        f = 1.0 / (2 * math.pi)
        J01, J1oz1 = J["1"]
        J02, J1oz2 = J["2"]
        k1, k2, b1, b2 = self.k1, self.k2, self.b1, self.b2
        # J1(kr) = k r J1oz
        A11 = -f * (b2 * k2 * k2 * J1oz2 - b1 * k1 * k1 * J1oz1) * r * md
        A12 = -f * (J01 - J02)
        A22 = f * (k1 * k1 * J1oz1 / b1 - k2 * k2 * J1oz2 / b2) * r * nd

        def AT(k, J0, J1oz):
            return -(k * k * f) * (J0 - J1oz) * nd * md - (k * k * f) * J1oz * (nm - nd * md)
        A21 = -(AT(k1, J01, J1oz1) - AT(k2, J02, J1oz2))
        return A11 / self.c1, A12 / self.c1, A21 / self.c2, A22 / self.c2

    def diag_B(self, nx, nz, kappa):
        """Diagonal limits of B = K - A log r for the periodic parts (scaled rows)."""
        g1, g2 = self.g1, self.g2
        f = 1.0 / (4 * math.pi)
        B11 = (self.b2 - self.b1) * kappa * f
        B22 = (1 / self.b1 - 1 / self.b2) * kappa * f
        B12 = (g1.spat_reg0 + g1.spec0) - (g2.spat_reg0 + g2.spec0)

        def BT(g, k):
            return (k * k / (8 * math.pi)) * (-EULER - 2 * math.log(g.E)) + (g.E ** 2 / (2 * math.pi)) * g.sum_q2 \
                - (nx * nx * g.spec_xx0 + nz * nz * g.spec_zz0)
        dBT = BT(g1, self.k1) - BT(g2, self.k2) - (g1.E ** 2 - g2.E ** 2) / (2 * math.pi)
        one = np.ones_like(kappa)
        return B11 / self.c1, B12 * one / self.c1, -dBT / self.c2, B22 / self.c2
