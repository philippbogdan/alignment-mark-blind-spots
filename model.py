"""The measured quantities and the indistinguishability rule of MODEL.md, and the reader of pair files."""
import math
import numpy as np

P = 2000.0
CHANNELS = [(532.0, "TE"), (532.0, "TM"), (633.0, "TE"), (633.0, "TM"),
            (775.0, "TE"), (775.0, "TM"), (850.0, "TE"), (850.0, "TM")]
NAMES = ["%d%s" % (lam, pol) for lam, pol in CHANNELS]
NOMINAL = np.array([1000.0, 160.0, 800.0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


def read(path):
    """A pair file as a dict: pair, cap, channels (names), sigma (position, visibility, strength;
    0.0 = must be equal, None = not read), kind, mark1, mark2, placement1, placement2, separation,
    and data (per channel: both marks' reading, visibility and strength as bie/ returns them)."""
    w = {"data": {}}
    for ln in open(path):
        ln = ln.split("#")[0].strip()
        if not ln:
            continue
        k, v = ln.split(None, 1)
        if k in ("data", "other"):  # per channel: X1 X2 V1 V2 S1 S2 (other: a channel outside the claim)
            name, *vals = v.split()
            w["data"][name] = [float(x) for x in vals]
            continue
        if k in w:
            raise ValueError("%s: field %r given twice" % (path, k))
        if k in ("mark1", "mark2"):
            w[k] = np.array([float(x) for x in v.split()])
            assert len(w[k]) == 12, k
        elif k == "channels":
            w[k] = v.split()
        elif k == "sigma":
            w[k] = [None if x == "-" else float(x) for x in v.split()]
        elif k in ("cap", "placement1", "placement2", "separation"):
            w[k] = float(v)
        else:
            w[k] = v
    return w


def mirror(p):
    """The mirror image x -> -x of a mark, as parameters."""
    w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = [float(v) for v in p]
    return np.array([w, h + b, t - b, cR, cL, aR, aL, -b, pR, pL, qR, qL])


def wrap(d):
    """Shortest signed distance modulo 1000 nm (readings are defined modulo half a pitch)."""
    return (np.asarray(d) + 500.0) % 1000.0 - 500.0


def measured(Ap, Am, ep, em, s0):
    """Shift of the position reading from the placement (nm), fringe visibility, signal strength
    and order imbalance from the +1 and -1 order amplitudes and efficiencies; s0 is e+ + e- of the
    nominal mark in the same channels."""
    shift = -(P / (4 * math.pi)) * np.angle(Ap / Am)
    V = 2 * np.sqrt(ep * em) / (ep + em)
    S = (ep + em) / s0
    B = (ep - em) / (ep + em)
    return shift, V, S, B


def compare(pair, X, V, S):
    """How far apart the two marks' data are on the pair's own channels.  X, V, S: arrays
    (2 marks, 8 channels).  Returns (d2, equal): d2 sums the squared differences of the noisy
    quantities over their variances (the two marks are indistinguishable when d2 <= 1); equal
    is the largest difference of each quantity that must agree exactly."""
    idx = [NAMES.index(c) for c in pair["channels"]]
    d2, equal = 0.0, {}
    for name, q, sig in zip(("position", "visibility", "strength"), (X, V, S), pair["sigma"]):
        if sig is None:
            continue
        diff = q[0, idx] - q[1, idx]
        if name == "position":
            diff = wrap(diff)
        if sig == 0.0:
            equal[name] = float(np.abs(diff).max())
        else:
            d2 += float(np.sum(diff ** 2)) / sig ** 2
    return d2, equal
