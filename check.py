"""Recompute a pair from its plain-text file and compare with its record.

    python check.py pairs/01-r80-all8.txt              boundary-integral solver (numpy, scipy)
    python check.py pairs/01-r80-all8.txt --fem        finite elements, order 8 (needs ngsolve)
    python check.py pairs/01-r80-all8.txt --fem 10     finite elements, order 10 (about 6 GB)
    python check.py pairs/01-r80-all8.txt --band 400 1100 10    readings across a band (csv);
                                                                with --fem, by finite elements

Both marks are solved in all eight channels; the channels the pair claims are marked with *.
"""
import json
import os
import sys

import numpy as np

import model as PR

HERE = os.path.dirname(os.path.abspath(__file__))


def admissible(pair):
    from fem.geometry import admissibility
    return [admissibility(pair[m], pair["cap"]) for m in ("mark1", "mark2")]


def solve(marks, channels, method, order=8):
    """A+, A-, e+, e- for the given marks and for the nominal mark (last row)."""
    marks = np.vstack([np.atleast_2d(marks), PR.NOMINAL])
    if method == "bie":
        from bie.solver import solve as bie_solve
        return bie_solve(marks, channels)
    from fem.solver import amplitudes
    return amplitudes(marks, channels, order=order)


def band(pair, lo, hi, step, order=None):
    """Shift and order efficiencies of mark1 (and of mark2 unless it is the mirror image, whose
    shift is the opposite and whose two orders are swapped) across a band, by boundary integrals
    or, given an order, by finite elements."""
    marks = np.array([pair["mark1"]] + ([] if pair["kind"] == "mirror" else [pair["mark2"]]))
    channels = [(lam, pol) for lam in np.arange(lo, hi + 1e-9, step) for pol in ("TE", "TM")]
    print("wavelength_nm,polarisation," + ",".join("shift%d_nm,eplus%d,eminus%d" % (k, k, k) for k in range(1, len(marks) + 1)))

    def row(lam, pol, sol, j):
        shift = -(PR.P / (4 * np.pi)) * np.angle(sol["Ap"][:, j] / sol["Am"][:, j])
        print("%g,%s," % (lam, pol) + ",".join("%.12f,%.12e,%.12e" % (shift[k], sol["ep"][k, j], sol["em"][k, j])
                                               for k in range(len(marks))), flush=True)

    if order is None:
        from bie.solver import solve as bie_solve
        for lam, pol in channels:
            row(lam, pol, bie_solve(marks, [(lam, pol)]), 0)
    else:
        from fem.solver import amplitudes
        sol = amplitudes(marks, channels, order=order)
        for j, (lam, pol) in enumerate(channels):
            row(lam, pol, sol, j)


def main(argv):
    if not argv or argv[0].startswith("-"):
        sys.exit(__doc__)
    path = argv[0]
    pair = PR.read(path)
    method, order = "bie", None
    if "--fem" in argv:
        i = argv.index("--fem")
        method, order = "fem", int(argv[i + 1]) if len(argv) > i + 1 and argv[i + 1].isdigit() else 8
    if "--band" in argv:
        i = argv.index("--band")
        return band(pair, *[float(v) for v in argv[i + 1:i + 4]], order=order)
    bad = admissible(pair)
    print("pair %s: cap %g nm, %s, placements %.10g and %.10g nm (%.6f nm apart)" % (
        pair["pair"], pair["cap"], pair["kind"], pair["placement1"], pair["placement2"], pair["separation"]))
    print("both marks within the cap: %s" % ("yes" if not (bad[0] or bad[1]) else "NO %s" % bad))
    if pair["kind"] == "mirror":
        print("mark2 is the mirror image of mark1: %s" % (
            "yes" if np.allclose(PR.mirror(pair["mark1"]), pair["mark2"], rtol=0, atol=1e-9) else "NO"))
    sol = solve([pair["mark1"], pair["mark2"]], PR.CHANNELS, method, order)
    s0 = sol["ep"][2] + sol["em"][2]
    shift, V, S, B = PR.measured(sol["Ap"][:2], sol["Am"][:2], sol["ep"][:2], sol["em"][:2], s0)
    X = shift + np.array([[pair["placement1"]], [pair["placement2"]]])
    print("\n%-8s %16s %16s %15s %15s %13s %13s %11s %11s" % (
        "channel", "reading 1 (nm)", "reading 2 (nm)", "visibility 1", "visibility 2", "strength 1", "strength 2",
        "imbalance 1", "imbalance 2"))
    for j, name in enumerate(PR.NAMES):
        print("%-8s %16.10f %16.10f %15.12f %15.12f %13.10f %13.10f %11.2e %11.2e" % (
            name + ("*" if name in pair["channels"] else " "), X[0, j], X[1, j], V[0, j], V[1, j], S[0, j], S[1, j],
            B[0, j], B[1, j]))
    d2, equal = PR.compare(pair, X, V, S)
    print()
    for name, val in equal.items():
        print("largest difference in %s over the %d claimed channels: %.1e" % (name, len(pair["channels"]), val))
    if any(s for s in pair["sigma"] if s):
        print("d^2 over the claimed channels with noise %s: %.4f (indistinguishable when at most 1)" % (
            " ".join("-" if s is None else "%g" % s for s in pair["sigma"]), d2))
    rec_path = os.path.join(HERE, "records", os.path.basename(path).replace(".txt", ".json"))
    if os.path.exists(rec_path):
        rec = json.load(open(rec_path))
        key = "bie" if method == "bie" else "fem%d" % order
        if key in rec:
            r = {k: np.array(rec[key][k], dtype=float) for k in ("shift", "visibility", "strength")}   # nan: not recorded
            dX = float(np.nanmax(np.abs(PR.wrap(shift - r["shift"]))))
            dV = float(np.nanmax(np.abs(V - r["visibility"])))
            dS = float(np.nanmax(np.abs(S - r["strength"])))
            print("against the record (%s): readings within %.1e nm, visibility %.1e, strength %.1e" % (key, dX, dV, dS))


if __name__ == "__main__":
    main(sys.argv[1:])
