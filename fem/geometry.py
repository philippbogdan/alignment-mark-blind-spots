"""Mark geometry exactly as fixed in MODEL.md.

Parameter vector (all lengths in nm), order used everywhere in this repository:

    0 w   line width            [950, 1050]
    1 h   line height (left top corner height)   [144, 176]
    2 t   oxide thickness (oxide plane at z = h + t)   [760, 840]
    3 cL  floor tilt, left foot height   [0, r]
    4 cR  floor tilt, right foot height  [0, r]
    5 aL  sidewall slope, left top corner inward shift  [0, r]
    6 aR  sidewall slope, right top corner inward shift [0, r]
    7 b   top tilt, right top corner at h + b   [-r, r]
    8 pL  rounding radius, left top corner   [0, r]
    9 pR  rounding radius, right top corner  [0, r]
   10 qL  chamfer, left top corner   [0, r]
   11 qR  chamfer, right top corner  [0, r]

The mark is placed at x0 = 0: feet at (-w/2, cL) and (w/2, cR).  Silicon fills everything
below the profile; the profile of one period is the line's upper boundary from the left foot
to the right foot plus the straight trench floor from the right foot (w/2, cR) to the next
line's left foot (P - w/2, cL).  At every height z the silicon in one period is a single
interval [xl(z), xr(z)] (or the whole period below both feet).
"""
import math
import numpy as np

P = 2000.0
NAMES = ["w", "h", "t", "cL", "cR", "aL", "aR", "b", "pL", "pR", "qL", "qR"]
NOMINAL = np.array([1000.0, 160.0, 800.0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
SYM_BOUNDS = {"w": (950.0, 1050.0), "h": (144.0, 176.0), "t": (760.0, 840.0)}


def mirror(p):
    """Mirror image x -> -x of a mark, as a parameter vector (placement maps x0 -> -x0)."""
    w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = [float(v) for v in p]
    return np.array([w, h + b, t - b, cR, cL, aR, aL, -b, pR, pL, qR, qL])


# ----------------------------------------------------------------------------------------
# primitives of the line's upper boundary
# ----------------------------------------------------------------------------------------
class Seg:
    """Straight segment from (x1, z1) to (x2, z2)."""

    def __init__(self, x1, z1, x2, z2):
        self.x1, self.z1, self.x2, self.z2 = x1, z1, x2, z2

    def zrange(self):
        return min(self.z1, self.z2), max(self.z1, self.z2)

    def x_at(self, z):
        if self.z2 == self.z1:
            raise ValueError("horizontal segment has no x(z)")
        return self.x1 + (self.x2 - self.x1) * (z - self.z1) / (self.z2 - self.z1)

    def x_avg(self, za, zb):
        return self.x_at(0.5 * (za + zb))


class Arc:
    """Arc piece of the circle (cx, cz, rad) on one side: x = cx + sgn*sqrt(rad^2-(z-cz)^2),
    for z in [zlo, zhi] (within the upper half of the circle)."""

    def __init__(self, cx, cz, rad, sgn, zlo, zhi):
        self.cx, self.cz, self.rad, self.sgn, self.zlo, self.zhi = cx, cz, rad, sgn, zlo, zhi

    def zrange(self):
        return self.zlo, self.zhi

    def x_at(self, z):
        u = (z - self.cz) / self.rad
        u = min(1.0, max(-1.0, u))
        return self.cx + self.sgn * self.rad * math.sqrt(1.0 - u * u)

    def x_avg(self, za, zb):
        r = self.rad

        def F(z):
            u = (z - self.cz) / r
            u = min(1.0, max(-1.0, u))
            return 0.5 * r * r * (u * math.sqrt(1.0 - u * u) + math.asin(u))

        return self.cx + self.sgn * (F(zb) - F(za)) / (zb - za)


def _unit(vx, vz):
    n = math.hypot(vx, vz)
    return vx / n, vz / n, n


class Profile:
    """Boundary of one mark, built from the parameter vector."""

    def __init__(self, p):
        p = [float(v) for v in p]
        self.p = p
        w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = p
        self.w, self.h, self.t, self.cL, self.cR = w, h, t, cL, cR
        FL = (-w / 2, cL)
        FR = (w / 2, cR)
        TL = (-w / 2 + aL, h)
        TR = (w / 2 - aR, h + b)
        self.FL, self.FR, self.TL, self.TR = FL, FR, TL, TR
        # unit vectors along the edges away from each top corner
        sLx, sLz, lenL = _unit(FL[0] - TL[0], FL[1] - TL[1])      # TL -> FL (down left sidewall)
        sRx, sRz, lenR = _unit(FR[0] - TR[0], FR[1] - TR[1])      # TR -> FR (down right sidewall)
        tx, tz, lenT = _unit(TR[0] - TL[0], TR[1] - TL[1])        # TL -> TR (along top)
        self.lenL, self.lenR, self.lenT = lenL, lenR, lenT
        # interior angles at the corners
        thL = math.acos(max(-1.0, min(1.0, sLx * tx + sLz * tz)))
        thR = math.acos(max(-1.0, min(1.0, sRx * (-tx) + sRz * (-tz))))
        self.thL, self.thR = thL, thR
        # distances of the tangent / cut points from the corner along each edge
        if pL > 0 and qL > 0:
            self.dL = float("nan")
        elif pL > 0:
            self.dL = pL / math.tan(thL / 2)
        else:
            self.dL = qL
        if pR > 0 and qR > 0:
            self.dR = float("nan")
        elif pR > 0:
            self.dR = pR / math.tan(thR / 2)
        else:
            self.dR = qR
        prims = []
        # left corner
        if self.dL > 0:
            A1 = (TL[0] + self.dL * sLx, TL[1] + self.dL * sLz)   # on sidewall
            A2 = (TL[0] + self.dL * tx, TL[1] + self.dL * tz)     # on top
        else:
            A1 = A2 = TL
        if self.dR > 0:
            B1 = (TR[0] + self.dR * sRx, TR[1] + self.dR * sRz)   # on right sidewall
            B2 = (TR[0] - self.dR * tx, TR[1] - self.dR * tz)     # on top
        else:
            B1 = B2 = TR
        self.A1, self.A2, self.B1, self.B2 = A1, A2, B1, B2
        prims.append(("seg", FL, A1))
        if self.dL > 0:
            if pL > 0:
                prims.append(("arc", self._arc_center(TL, (sLx, sLz), (tx, tz), pL, thL), pL, A1, A2))
            else:
                prims.append(("seg", A1, A2))
        prims.append(("seg", A2, B2))
        if self.dR > 0:
            if pR > 0:
                prims.append(("arc", self._arc_center(TR, (sRx, sRz), (-tx, -tz), pR, thR), pR, B2, B1))
            else:
                prims.append(("seg", B2, B1))
        prims.append(("seg", B1, FR))
        self.prims = prims
        self._build_monotone()

    @staticmethod
    def _arc_center(C, e1, e2, rad, th):
        bx, bz = e1[0] + e2[0], e1[1] + e2[1]
        bn = math.hypot(bx, bz)
        dist = rad / math.sin(th / 2)
        return (C[0] + dist * bx / bn, C[1] + dist * bz / bn)

    def _build_monotone(self):
        """Split the upper boundary into an ascending part (x increasing with z) and a
        descending part (x decreasing with z), as lists of (zlo, zhi, primitive)."""
        asc, desc = [], []
        apex_z = -1e30
        pieces = []   # (z_start, z_end, primitive) in traversal order, each monotone in z
        for pr in self.prims:
            if pr[0] == "seg":
                (x1, z1), (x2, z2) = pr[1], pr[2]
                if abs(x2 - x1) < 1e-15 and abs(z2 - z1) < 1e-15:
                    continue
                pieces.append((z1, z2, Seg(x1, z1, x2, z2)))
            else:
                (cx, cz), rad, S, E = pr[1], pr[2], pr[3], pr[4]
                # traversal from S to E goes clockwise over the top of the circle
                phS = math.atan2(S[1] - cz, S[0] - cx)
                phE = math.atan2(E[1] - cz, E[0] - cx)
                # clockwise traversal: angle decreases from phS to phE
                if phS > math.pi / 2 > phE:
                    top = cz + rad
                    pieces.append((S[1], top, Arc(cx, cz, rad, -1.0, S[1], top)))
                    pieces.append((top, E[1], Arc(cx, cz, rad, +1.0, E[1], top)))
                else:
                    sgn = -1.0 if (phS + phE) / 2 > math.pi / 2 else +1.0
                    pieces.append((S[1], E[1], Arc(cx, cz, rad, sgn, min(S[1], E[1]), max(S[1], E[1]))))
        for (za, zb, prim) in pieces:
            apex_z = max(apex_z, za, zb)
        self.apex = apex_z
        # ascending: pieces before the apex with zb > za ; descending: after apex with zb < za
        reached = False
        for (za, zb, prim) in pieces:
            if not reached:
                if zb > za:
                    asc.append((za, zb, prim))
                if max(za, zb) >= apex_z:
                    reached = True
                    if zb < za:  # starts at apex and descends
                        desc.append((zb, za, prim))
            else:
                if zb < za:
                    desc.append((zb, za, prim))
                elif zb > za + 1e-12:
                    raise ValueError("boundary not unimodal")
        self.asc, self.desc = asc, desc

    # ------------------------------------------------------------------------------
    def breakpoints(self):
        zs = {self.cL, self.cR, self.apex}
        for (za, zb, _) in self.asc + self.desc:
            zs.add(za)
            zs.add(zb)
        zmin = min(self.cL, self.cR)
        return sorted(z for z in zs if z >= zmin - 1e-12)

    def _find(self, lst, z):
        for (za, zb, prim) in lst:
            if za - 1e-12 <= z <= zb + 1e-12:
                return prim
        raise ValueError("height %r outside boundary part" % z)

    def interval_avg(self, za, zb):
        """Silicon interval [xl, xr] averaged over the layer za < z < zb (exact area).
        Requires that no breakpoint lies strictly inside (za, zb).  Returns None if the
        layer is entirely above the mark, 'full' if entirely silicon."""
        w, cL, cR = self.w, self.cL, self.cR
        zm = 0.5 * (za + zb)
        if zm >= self.apex:
            return None
        if zm < min(cL, cR):
            return "full"
        # left end
        if zm > cL:
            xl = self._find(self.asc, zm).x_avg(za, zb)
        else:
            # below the left foot: floor crossing in the trench to the left (floor rises to cL)
            # floor from (-P + w/2, cR) to (-w/2, cL), here cR < z < cL
            xl = (-P + w / 2) + (zm - cR) * (P - w) / (cL - cR)
        if zm > cR:
            xr = self._find(self.desc, zm).x_avg(za, zb)
        else:
            # below the right foot: floor from (w/2, cR) to (P - w/2, cL) with cL < z < cR
            xr = (w / 2) + (cR - zm) * (P - w) / (cR - cL)
        return xl, xr

    def _x_point(self, z, side):
        w, cL, cR = self.w, self.cL, self.cR
        if z >= self.apex or z < min(cL, cR):
            return None
        if side == 0:
            if z > cL:
                return self._find(self.asc, z).x_at(z)
            return (-P + w / 2) + (z - cR) * (P - w) / (cL - cR)
        if z > cR:
            return self._find(self.desc, z).x_at(z)
        return (w / 2) + (cR - z) * (P - w) / (cR - cL)


# ----------------------------------------------------------------------------------------
def admissibility(p, r, tol=0.0):
    """Return list of violated rules (empty list = admissible) for cap r."""
    p = [float(v) for v in p]
    w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = p
    bad = []
    for name, v in zip(NAMES[:3], p[:3]):
        lo, hi = SYM_BOUNDS[name]
        if not (lo - tol <= v <= hi + tol):
            bad.append("%s out of range" % name)
    for name, v in zip(NAMES[3:], p[3:]):
        lo = -r if name == "b" else 0.0
        if not (lo - tol <= v <= r + tol):
            bad.append("%s out of [%g, %g]" % (name, lo, r))
    if pL > 0 and qL > 0:
        bad.append("left corner has both rounding and chamfer")
    if pR > 0 and qR > 0:
        bad.append("right corner has both rounding and chamfer")
    if not (h > cL + tol * 0):
        bad.append("left sidewall does not rise")
    if not (h + b > cR):
        bad.append("right sidewall does not rise")
    if bad:
        return bad
    pr = Profile(p)
    if pr.dL > pr.lenL + tol:
        bad.append("left corner point beyond left sidewall")
    if pr.dR > pr.lenR + tol:
        bad.append("right corner point beyond right sidewall")
    if pr.dL + pr.dR > pr.lenT + tol:
        bad.append("corner modifications overlap on the top")
    return bad


def is_admissible(p, r, tol=0.0):
    return len(admissibility(p, r, tol)) == 0


# ----------------------------------------------------------------------------------------
def layers(p, dz_max=2.0, dx_max=8.0, nmin=1):
    """Staircase of the grating region of mark p.

    Returns (profile, thick, xl, xr, rep): from the bottom (z = min(cL, cR)) to the apex, each
    entry is a block of `rep` identical sublayers of thickness `thick` with silicon on the
    exact-area averaged interval [xl, xr].  Each interval between consecutive breakpoints is
    split uniformly into max(nmin, ceil(dz/dz_max), ceil(dx/dx_max)) sublayers, dx being the
    variation of the interval ends over it; sublayers with identical cross-sections (vertical
    walls) are merged into one block."""
    pr = Profile(p)
    zs = pr.breakpoints()
    thick, xl, xr, rep = [], [], [], []
    for z0, z1 in zip(zs[:-1], zs[1:]):
        if z1 - z0 <= 1e-12:
            continue
        dxv = 0.0
        for side in (0, 1):
            a = pr._x_point(z0 + 1e-9 * (z1 - z0), side)
            bq = pr._x_point(z1 - 1e-9 * (z1 - z0), side)
            if a is not None and bq is not None:
                dxv = max(dxv, abs(bq - a))
        n = max(nmin, int(math.ceil((z1 - z0) / dz_max)), int(math.ceil(dxv / dx_max)))
        edges = np.linspace(z0, z1, n + 1)
        ivs = [pr.interval_avg(za, zb) for za, zb in zip(edges[:-1], edges[1:])]
        if dxv < 1e-9 and all(iv is not None and iv != "full" for iv in ivs):
            # identical cross-sections: one block
            thick.append((z1 - z0) / n)
            xl.append(ivs[0][0])
            xr.append(ivs[0][1])
            rep.append(n)
            continue
        for (za, zb), iv in zip(zip(edges[:-1], edges[1:]), ivs):
            if iv is None:
                continue
            thick.append(zb - za)
            if iv == "full":
                xl.append(-P / 2)
                xr.append(P / 2)
            else:
                xl.append(iv[0])
                xr.append(iv[1])
            rep.append(1)
    return pr, np.array(thick), np.array(xl), np.array(xr), np.array(rep, dtype=int)
