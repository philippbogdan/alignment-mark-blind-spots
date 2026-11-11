"""Profile of one period of the alignment mark (own construction from MODEL.md).

params = [w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR], placement x0 = 0, period P = 2000.

The interface Gamma (Si below, SiO2 above) is traversed left to right over one period:
left foot (-w/2, cL) -> left sidewall -> [left corner chamfer/arc] -> top -> [right corner
chamfer/arc] -> right sidewall -> right foot (w/2, cR) -> floor -> next left foot (P - w/2, cL).
Unit normal n = left normal of the traversal direction (points from Si into SiO2).
Every junction between consecutive segments is a breakpoint (corner, or C1 tangent point of an
arc); each is resolved by RCIP.
"""
import math
import numpy as np

P = 2000.0


class Line:
    kind = "line"

    def __init__(self, p0, p1):
        self.p0 = np.asarray(p0, float)
        self.p1 = np.asarray(p1, float)
        v = self.p1 - self.p0
        self.length = float(np.hypot(*v))
        self.tau = v / self.length

    def at(self, s):
        """absolute position, tangent, curvature at arclength s from start (array)."""
        s = np.asarray(s, float)
        x = self.p0[0] + s * self.tau[0]
        z = self.p0[1] + s * self.tau[1]
        one = np.ones_like(s)
        return x, z, self.tau[0] * one, self.tau[1] * one, 0.0 * one

    def rel_start(self, sig):
        """offset from the start point at arclength sig from the start."""
        sig = np.asarray(sig, float)
        one = np.ones_like(sig)
        return sig * self.tau[0], sig * self.tau[1], self.tau[0] * one, self.tau[1] * one, 0 * one

    def rel_end(self, sig):
        """offset from the end point at arclength sig before the end."""
        sig = np.asarray(sig, float)
        one = np.ones_like(sig)
        return -sig * self.tau[0], -sig * self.tau[1], self.tau[0] * one, self.tau[1] * one, 0 * one


class Arc:
    kind = "arc"

    def __init__(self, center, radius, phi0, phi1):
        """circular arc from angle phi0 to phi1 (traversed in that direction, |phi1-phi0| < pi)."""
        self.c = np.asarray(center, float)
        self.R = float(radius)
        self.phi0, self.phi1 = float(phi0), float(phi1)
        self.sgn = 1.0 if phi1 > phi0 else -1.0
        self.length = self.R * abs(phi1 - phi0)
        self.kappa = self.sgn / self.R
        self.p0 = self.c + self.R * np.array([math.cos(phi0), math.sin(phi0)])
        self.p1 = self.c + self.R * np.array([math.cos(phi1), math.sin(phi1)])

    def _tan(self, phi):
        return -self.sgn * np.sin(phi), self.sgn * np.cos(phi)

    def at(self, s):
        s = np.asarray(s, float)
        phi = self.phi0 + self.sgn * s / self.R
        tx, tz = self._tan(phi)
        return (self.c[0] + self.R * np.cos(phi), self.c[1] + self.R * np.sin(phi), tx, tz,
                self.kappa * np.ones_like(s))

    def _rel(self, phiA, dphi):
        # p(phiA + dphi) - p(phiA), stable for small dphi
        sh = np.sin(dphi / 2)
        m = phiA + dphi / 2
        return -2 * self.R * sh * np.sin(m), 2 * self.R * sh * np.cos(m)

    def rel_start(self, sig):
        sig = np.asarray(sig, float)
        dphi = self.sgn * sig / self.R
        dx, dz = self._rel(self.phi0, dphi)
        tx, tz = self._tan(self.phi0 + dphi)
        return dx, dz, tx, tz, self.kappa * np.ones_like(sig)

    def rel_end(self, sig):
        sig = np.asarray(sig, float)
        dphi = -self.sgn * sig / self.R
        dx, dz = self._rel(self.phi1, dphi)
        tx, tz = self._tan(self.phi1 + dphi)
        return dx, dz, tx, tz, self.kappa * np.ones_like(sig)


def _unit(v):
    v = np.asarray(v, float)
    return v / np.hypot(*v)


def _corner(C, u_in, u_out, p, q):
    """Modify the corner C, reached along -u_in... Here u_in points from C back along the
    incoming edge, u_out from C along the outgoing edge.  Returns (point_on_incoming_edge,
    list_of_segments_between, point_on_outgoing_edge, distance_used_on_each_edge)."""
    if p > 0 and q > 0:
        raise ValueError("corner carries both rounding and chamfer")
    if q > 0:
        A = C + q * u_in
        B = C + q * u_out
        return A, [Line(A, B)], B, q
    if p > 0:
        cth = float(np.clip(np.dot(u_in, u_out), -1, 1))
        th = math.acos(cth)                       # interior angle (Si side)
        d = p / math.tan(th / 2)
        A = C + d * u_in
        B = C + d * u_out
        bis = _unit(u_in + u_out)
        cen = C + (p / math.sin(th / 2)) * bis
        phiA = math.atan2(A[1] - cen[1], A[0] - cen[0])
        phiB = math.atan2(B[1] - cen[1], B[0] - cen[0])
        dphi = (phiB - phiA + math.pi) % (2 * math.pi) - math.pi
        return A, [Arc(cen, p, phiA, phiA + dphi)], B, d
    return C, [], C, 0.0


class Profile:
    def __init__(self, params, P_=P, x0=0.0):
        w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = [float(v) for v in params]
        self.params = np.array(params, float)
        self.P = P_
        self.H = h + t
        FL = np.array([x0 - w / 2, cL])
        FR = np.array([x0 + w / 2, cR])
        TL = np.array([x0 - w / 2 + aL, h])
        TR = np.array([x0 + w / 2 - aR, h + b])
        FLn = np.array([x0 + P_ - w / 2, cL])
        # left top corner: incoming edge = left sidewall, outgoing = top
        A1, mid1, B1, d1 = _corner(TL, _unit(FL - TL), _unit(TR - TL), pL, qL)
        A2, mid2, B2, d2 = _corner(TR, _unit(TL - TR), _unit(FR - TR), pR, qR)
        segs = [Line(FL, A1)] + mid1 + [Line(B1, A2)] + mid2 + [Line(B2, FR), Line(FR, FLn)]
        self.segs = segs
        # admissibility (the parts we rely on)
        side_l = np.hypot(*(TL - FL))
        side_r = np.hypot(*(TR - FR))
        top = np.hypot(*(TR - TL))
        ok = (h > cL) and (h + b > cR) and d1 < side_l and d2 < side_r and d1 + d2 < top \
            and all(s.length > 0 for s in segs)
        if not ok:
            raise ValueError("inadmissible profile %s" % (params,))
        zs = []
        for s in segs:
            ss = np.linspace(0, s.length, 65)
            zs.append(s.at(ss)[1])
        zs = np.concatenate(zs)
        self.zmin, self.zmax = float(zs.min()), float(zs.max())
        if self.zmax >= self.H:
            raise ValueError("profile reaches the oxide surface")

    def corner_angles(self):
        """interior angle on the oxide side at each breakpoint j (between segs[j-1] and segs[j])."""
        out = []
        n = len(self.segs)
        for j in range(n):
            a, b = self.segs[j - 1], self.segs[j]
            ta = np.array(a.at(a.length)[2:4]).ravel()
            tb = np.array(b.at(0.0)[2:4]).ravel()
            turn = math.atan2(ta[0] * tb[1] - ta[1] * tb[0], ta @ tb)   # left turn > 0
            out.append(math.pi - turn)
        return np.array(out)


def mirror_params(p):
    """Mirror image x -> -x of a mark, expressed in the same parametrisation."""
    w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = p
    return np.array([w, h + b, t - b, cR, cL, aR, aL, -b, pR, pL, qR, qL], float)
