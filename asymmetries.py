"""What each kind of left-right asymmetry does to the eight readings, on the line of pair 01.

The line of pair 01 is first made symmetric (each left and right parameter replaced by their
mean, the top levelled).  Then one asymmetry at a time is put back, 0.5 nm between left and right:
the floor (left foot up, right foot down), the top (tilted), the walls (the left leaning more)
and the shoulders (the left rounder).  The last case is mark 1 of pair 01 itself, which mixes
the four.  Solver: bie/.  About ten minutes.

    python asymmetries.py > records/asymmetries.json
"""
import json
import sys

import numpy as np

import model as PR
from bie.solver import solve


def cases(mark, size=0.25):
    w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = mark
    line = np.array([w, h + b / 2, t - b / 2, (cL + cR) / 2, (cL + cR) / 2, (aL + aR) / 2, (aL + aR) / 2, 0.0,
                     (pL + pR) / 2, (pL + pR) / 2, 0.0, 0.0])

    def with_(floor=0.0, wall=0.0, shoulder=0.0, top=0.0):
        q = line.copy()
        q[3] += floor; q[4] -= floor
        q[5] += wall; q[6] -= wall
        q[8] += shoulder; q[9] -= shoulder
        q[7] = top; q[1] -= top / 2; q[2] += top / 2
        return q

    return [("symmetric", "the line of pair 01 with its left and right sides made equal", line),
            ("floor", "left foot 0.25 nm higher, right foot 0.25 nm lower", with_(floor=size)),
            ("top", "top tilted, its right end 0.5 nm higher than its left", with_(top=2 * size)),
            ("wall", "left wall leaning 0.25 nm more, right wall 0.25 nm less", with_(wall=size)),
            ("shoulder", "left shoulder radius 0.25 nm larger, right 0.25 nm smaller", with_(shoulder=size)),
            ("pair 01", "mark 1 of pair 01: floor %+.4f, top %+.4f, wall %+.4f, shoulder %+.4f nm between left and right" % (
                cL - cR, b, aL - aR, pL - pR), np.array(mark, float))]


if __name__ == "__main__":
    pair = PR.read("pairs/01-r80-all8-identical.txt")
    out = dict(channels=PR.NAMES, cases=[])
    for name, what, mark in cases(pair["mark1"]):
        sol = solve(mark, PR.CHANNELS)
        shift = -(PR.P / (4 * np.pi)) * np.angle(sol["Ap"][0] / sol["Am"][0])
        out["cases"].append(dict(name=name, what=what, mark=mark.tolist(), shift=shift.tolist()))
        print(name, np.round(shift, 6), file=sys.stderr, flush=True)
    print(json.dumps(out, indent=1))
