"""Draws the figures of the README from the files in records/ and pairs/.

    python figures/make.py        (matplotlib)
"""
import csv
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import model as PR                      # noqa: E402
from fem.geometry import Profile        # noqa: E402

plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"], "font.size": 10.5,
                     "svg.fonttype": "path", "svg.hashsalt": "figures", "axes.linewidth": 0.8, "axes.edgecolor": "#444444",
                     "xtick.color": "#444444", "ytick.color": "#444444", "text.color": "#1f2328",
                     "axes.labelcolor": "#1f2328"})
INK, MUTED, RULE = "#1f2328", "#6b7280", "#d0d7de"
SILICON, OXIDE = "#8d99a6", "#dbeafe"
COLOUR = {532: "#1a9641", 633: "#e6550d", 775: "#a50f15", 850: "#525252"}
HEAD = "01-r80-all8-identical"


def outline(mark, x0=0.0, n=200):
    """The line's upper boundary from the left foot to the right foot."""
    pr = Profile(mark)
    pts = []
    for prim in pr.prims:
        if prim[0] == "seg":
            pts += [prim[1], prim[2]]
        else:
            (cx, cz), rad, s, e = prim[1], prim[2], prim[3], prim[4]
            a0, a1 = math.atan2(s[1] - cz, s[0] - cx), math.atan2(e[1] - cz, e[0] - cx)
            da = (a1 - a0 + math.pi) % (2 * math.pi) - math.pi
            pts += [(cx + rad * math.cos(a0 + da * k / n), cz + rad * math.sin(a0 + da * k / n)) for k in range(n + 1)]
    pts = np.array(pts)
    pts[:, 0] += x0
    return pts


def ruler(height, x1, x2, left=46.0, right=72.0):
    """An 830 px wide figure (100 px to the inch) with one axes over all of it: x is position in nm, with x1
    `left` px from the left edge and x2 `right` px from the right edge; y is px from the bottom.  Returns the
    figure, the axes and px per nm."""
    ppn = (830.0 - left - right) / (x2 - x1)
    lo = x1 - left / ppn
    fig = plt.figure(figsize=(8.3, height / 100.0), dpi=100)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(lo, lo + 830.0 / ppn)
    ax.set_ylim(0.0, height)
    ax.axis("off")
    return fig, ax, ppn


def ruler_frame(ax, ppn, x1, x2, xr, names, count, top, ya):
    """What both reading figures share: the two marks as named vertical lines, the reading that `count`
    channels give both of them as a grey line, the distance from it to each mark, and the position axis at ya."""
    for x, name in ((x1, names[0]), (x2, names[1])):
        ax.plot([x, x], [ya, top - 8], color=INK, lw=1.6, solid_capstyle="butt", zorder=3)
        ax.text(x, top, name, ha="center", va="baseline", fontsize=11.5, fontweight="bold")
    ax.plot([xr, xr], [ya, top - 8], color=MUTED, lw=0.9, zorder=2)
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight"}
    ax.text(xr, top, "both marks give these %s readings" % words[count], ha="center", va="baseline",
            fontsize=10.5)
    yd = ya + 18
    arrow = dict(arrowstyle="<|-|>", color=INK, lw=1.3, shrinkA=0, shrinkB=0, mutation_scale=10)
    for a, b in ((x1, xr), (xr, x2)):
        ax.annotate("", (a, yd), (b, yd), arrowprops=arrow)
        ax.text(0.5 * (a + b), yd + 6, "%.4f nm" % (b - a), ha="center", va="baseline", fontsize=12.5,
                fontweight="bold")
    ticks = np.arange(0.0, max(x1, x2) + 0.05, 0.1)
    ax.plot([min(x1, 0.0) - 24 / ppn, max(x2, ticks[-1]) + 24 / ppn], [ya, ya], color="#444444", lw=0.8)
    for t in ticks:
        ax.plot([t, t], [ya, ya - 4], color="#444444", lw=0.8)
        ax.text(t, ya - 7, "0" if abs(t) < 1e-12 else "%.1f" % t, ha="center", va="top", fontsize=9.5, color=MUTED)
    ax.text(0.5 * (x1 + x2), ya - 26, "position (nm)", ha="center", va="top", fontsize=9.5, color=MUTED)


def ruler_reading(ax, ppn, x, y, channel, label_at=None):
    """One reading: a dot in the colour of its wavelength, filled for TE and open for TM.  Given label_at, the
    channel is named to the right of that position, and the name returned."""
    lam, pol = int(channel[:3]), channel[3:]
    c = COLOUR[lam]
    ax.scatter([x], [y], s=46, marker="o", facecolor=c if pol == "TE" else "white", edgecolor=c, linewidth=1.5,
               zorder=5)
    if label_at is None:
        return None
    return ax.text(label_at + 9.0 / ppn, y, "%d nm %s" % (lam, pol), ha="left", va="center", fontsize=10, color=c)


def figure_readings():
    """Pair 01: all eight channels read the mark and its mirror image at one position, 0.1998 nm from each."""
    name = "01-r80-all8-identical"
    pair = PR.read(os.path.join(ROOT, "pairs", name + ".txt"))
    rec = json.load(open(os.path.join(ROOT, "records", name + ".json")))
    x1, x2 = pair["placement1"], pair["placement2"]
    shift = np.array(rec["bie"]["shift"])
    reads = np.concatenate([x1 + shift[0], x2 + shift[1]])
    assert np.ptp(reads) < 1e-9, "the eight channels do not read one position"
    xr = float(np.mean(reads))
    chans = rec["channels"]
    # layout in px from the top: names, the channels every 19 px, the distances, the axis
    rows = [54.0 + 19.0 * k for k in range(len(chans))]
    axis = rows[-1] + 58.0
    H = axis + 52.0
    fig, ax, ppn = ruler(H, x1, x2)
    ruler_frame(ax, ppn, x1, x2, xr, ("the mark", "its mirror image"), len(chans), H - 24.0, H - axis)
    for j, ch in enumerate(chans):
        ruler_reading(ax, ppn, x1 + shift[0][j], H - rows[j], ch, label_at=xr)
    save(fig, "readings")


def spectrum_readings(lo=500, hi=1000):
    """shift1_nm of records/spectrum-01.csv from lo to hi nm, as {polarisation: (wavelengths, readings)}."""
    rows = list(csv.DictReader(open(os.path.join(ROOT, "records", "spectrum-01.csv"))))
    out = {}
    for pol in ("TE", "TM"):
        r = sorted((float(x["wavelength_nm"]), float(x["shift1_nm"])) for x in rows
                   if x["polarisation"] == pol and lo <= float(x["wavelength_nm"]) <= hi)
        out[pol] = (np.array([w for w, _ in r]), np.array([v for _, v in r]))
    return out


def finish_spectrum():
    """Mark 1 of pair 01 read at every wavelength from 500 to 1000 nm, twelve colours 45 nm apart picked out."""
    from matplotlib.transforms import blended_transform_factory
    lo, hi, step = 500, 1000, 45
    curves = spectrum_readings(lo, hi)
    lam = curves["TE"][0]
    assert np.array_equal(lam, curves["TM"][0]) and np.array_equal(lam, np.arange(lo, hi + 1)), "1 nm steps, TE and TM"
    picks = lam[(lam - lo) % step == 0]
    at = {pol: v[(lam - lo) % step == 0] for pol, (_, v) in curves.items()}
    every = np.concatenate([at["TE"], at["TM"]])
    rec = json.load(open(os.path.join(ROOT, "records", HEAD + ".json")))
    d0 = float(np.mean(rec["bie"]["shift"][0]))
    four = sorted({int(c[:-2]) for c in rec["channels"]})
    for w in four:  # the two curves meet there, at the reading of the first figure
        assert all(abs(curves[pol][1][lam == w][0] - d0) < 1e-9 for pol in ("TE", "TM")), w
    words = "zero one two three four five six seven eight nine ten eleven twelve".split()
    tm_ink = "#8250df"
    fig = plt.figure(figsize=(8.3, 3.4))
    ax = fig.add_axes([0.085, 0.155, 0.735, 0.81])
    for pol, c in (("TE", INK), ("TM", tm_ink)):
        ax.plot(lam, curves[pol][1], color=c, lw=1.1, zorder=3, solid_joinstyle="round")
    ax.scatter(picks, at["TM"], s=70, facecolor="white", edgecolor=tm_ink, linewidth=1.6, zorder=5, clip_on=False)
    ax.scatter(picks, at["TE"], s=46, facecolor=INK, edgecolor="white", linewidth=0.8, zorder=6, clip_on=False)
    # the four wavelengths of the first figure, a pointer under each
    for w in four:
        ax.plot([w, w], [d0 - 0.036, d0 - 0.027], color=COLOUR[w], lw=1.6, solid_capstyle="butt", zorder=4)
        ax.text(w, d0 - 0.041, "%d" % w, ha="center", va="top", fontsize=9, color=COLOUR[w])
    # the gap: from the mark to the readings
    xa = 0.5 * (picks[3] + picks[4])
    ax.annotate("", (xa, 0.0), (xa, every.mean()), arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=1.2, mutation_scale=11,
                                                                      shrinkA=0, shrinkB=0), zorder=2)
    ax.text(xa + 9, 0.5 * every.mean(), "%.1f nm" % every.mean(), ha="left", va="center", fontsize=15, fontweight="bold")
    # the spread of the readings
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    bx = 1.025
    ax.plot([bx - 0.009, bx, bx, bx - 0.009], [every.max(), every.max(), every.min(), every.min()], color=INK, lw=1.0,
            transform=tr, clip_on=False)
    ax.text(bx + 0.016, 0.5 * (every.max() + every.min()), "all %d readings\nwithin %.2f nm\nof each other"
            % (every.size, math.ceil(100 * (every.max() - every.min())) / 100), transform=tr, ha="left", va="center",
            fontsize=10.5, fontweight="bold",
            linespacing=1.3)
    ax.text(lo + 4, 0.007, "the mark", ha="left", va="bottom", fontsize=10.5, fontweight="bold")
    ax.text(hi, 0.262, "%s colours, %d nm apart, both polarisations" % (words[picks.size], np.diff(picks)[0]), ha="right",
            va="bottom", fontsize=10.5, color=INK)
    ax.text(915, 0.214, "TE", ha="center", va="bottom", fontsize=10.5, color=INK)
    ax.text(915, 0.175, "TM", ha="center", va="top", fontsize=10.5, color=tm_ink)
    ax.set_xlim(lo, hi)
    ax.set_ylim(0.0, 0.3)
    ax.set_yticks(np.arange(0.0, 0.31, 0.1))
    ax.set_yticklabels(["%g" % round(v, 1) for v in ax.get_yticks()])
    ax.set_xticks(np.arange(lo, hi + 1, 100))
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("position (nm)")
    ax.tick_params(length=3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_position(("outward", 10))
    ax.spines["bottom"].set_linewidth(1.8)
    ax.spines["bottom"].set_color(INK)
    save(fig, "spectrum")


def figure_small():
    """Pair 19: six channels read the line and its mirror image at one position; two tell them apart."""
    name = "19-r10-012467-identical"
    pair = PR.read(os.path.join(ROOT, "pairs", name + ".txt"))
    rec = json.load(open(os.path.join(ROOT, "records", name + ".json")))
    x1, x2 = pair["placement1"], pair["placement2"]
    shift = np.array(rec["bie"]["shift"])
    r1, r2 = x1 + shift[0], x2 + shift[1]
    chans = rec["channels"]
    same = [j for j in range(len(chans)) if rec["claimed"][j]]
    apart = [j for j in range(len(chans)) if not rec["claimed"][j]]
    reads = np.concatenate([r1[same], r2[same]])
    assert np.ptp(reads) < 1e-9, "the claimed channels do not read one position"
    xr = float(np.mean(reads))
    # layout in px from the top: as figure_readings, with the channels that tell the marks apart below the rest
    rows = [54.0 + 19.0 * k for k in range(len(same))]
    tags = rows[-1] + 38.0
    rows2 = [tags + 18.0 + 21.0 * k for k in range(len(apart))]
    axis = rows2[-1] + 62.0
    H = axis + 52.0
    fig, ax, ppn = ruler(H, x1, x2)
    ruler_frame(ax, ppn, x1, x2, xr, ("the line", "its mirror image"), len(same), H - 24.0, H - axis)
    for j, y in zip(same, rows):
        ruler_reading(ax, ppn, r1[j], H - y, chans[j], label_at=xr)
    # each mark's own reading in the channels that tell them apart, the two joined, and above them, on the side
    # of the common reading where they fall, whose readings they are
    for r, who in ((r1, "the line"), (r2, "its mirror image")):
        left = float(np.mean(r[apart])) < xr
        ax.text(xr + (-12.0 if left else 12.0) / ppn, H - tags, "readings of " + who, ha="right" if left else "left",
                va="baseline", fontsize=10)
    for j, y in zip(apart, rows2):
        ax.plot([r1[j], r2[j]], [H - y, H - y], color=COLOUR[int(chans[j][:3])], lw=1.5, zorder=4)
        ruler_reading(ax, ppn, min(r1[j], r2[j]), H - y, chans[j])
        t = ruler_reading(ax, ppn, max(r1[j], r2[j]), H - y, chans[j], label_at=max(r1[j], r2[j]))
        ax.annotate("%.2f nm apart" % abs(r1[j] - r2[j]), (1.0, 0.5), xycoords=t, xytext=(5, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=10, fontweight="bold")
    save(fig, "small")


def figure_mark():
    """The twelve parameters on an exaggerated line (not a mark of any pair, not to scale)."""
    m = np.array([760.0, 230.0, 250.0, 50.0, 26.0, 150.0, 95.0, 42.0, 120.0, 0.0, 0.0, 95.0])
    w, h, t, cL, cR, aL, aR, b, pL, pR, qL, qR = m
    pr = Profile(m)
    H = h + t
    X0, X1, Z0 = -640.0, 640.0, -150.0
    fig = plt.figure(figsize=(8.3, 4.28))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    slope = (cL - cR) / (2000 - w)
    floorL = [(X0, cL - slope * (-w / 2 - X0)), (-w / 2, cL)]
    floorR = [(w / 2, cR), (X1, cR + slope * (X1 - w / 2))]
    poly = np.vstack([[(X0, Z0)], floorL, outline(m), floorR, [(X1, Z0)]])
    ax.add_patch(Polygon([(X0, Z0), (X1, Z0), (X1, H), (X0, H)], closed=True, facecolor=OXIDE, edgecolor="none"))
    ax.add_patch(Polygon(poly, closed=True, facecolor=SILICON, edgecolor="none"))
    ax.plot([X0, X1], [H, H], color="#93b4dd", lw=1.0)
    ax.plot([X0, X1], [0, 0], color="white", lw=0.8, ls=(0, (4, 3)))
    ax.text(X0 + 12, -8, "z = 0", color="white", fontsize=9.5, va="top")
    ax.text(X0 + 12, H - 10, "z = h + t", color="#4b6a94", fontsize=9.5, va="top")
    ax.text(X0 + 12, Z0 + 22, "silicon", color="white", fontsize=10.5)
    ax.text(X0 + 12, 0.5 * (h + H), "oxide", color="#4b6a94", fontsize=10.5, va="center")
    ax.text(X0 + 12, H + 16, "air", color=MUTED, fontsize=10.5)
    ax.text(X1 - 12, H + 16, "not to scale: every deformation is drawn many times too large", ha="right", color=MUTED, fontsize=10.5)
    mono = dict(fontsize=12, family="Menlo")
    dash = dict(lw=0.8, ls=(0, (3, 3)))

    def dim(p0, p1, label, dx=0.0, dy=0.0, ha="center", va="center", colour=INK):
        ax.annotate("", p0, p1, arrowprops=dict(arrowstyle="<->", color=colour, lw=1.0, shrinkA=0, shrinkB=0))
        ax.text(0.5 * (p0[0] + p1[0]) + dx, 0.5 * (p0[1] + p1[1]) + dy, label, ha=ha, va=va, color=colour, **mono)

    FL, FR, TL, TR = pr.FL, pr.FR, pr.TL, pr.TR
    # width and placement
    ax.plot([FL[0], FL[0]], [-112, cL], color="white", **dash)
    ax.plot([FR[0], FR[0]], [-112, cR], color="white", **dash)
    dim((FL[0], -95), (FR[0], -95), "w", dx=110, dy=-20, colour="white")
    ax.plot([0, 0], [-112, -78], color="white", lw=1.6)
    ax.text(0, -70, "x0", ha="center", va="bottom", color="white", **mono)
    # floor tilt
    dim((FL[0] - 70, 0), (FL[0] - 70, cL), "cL", dx=-14, ha="right")
    dim((FR[0] + 70, 0), (FR[0] + 70, cR), "cR", dx=14, ha="left")
    ax.plot([FL[0] - 90, FL[0]], [cL, cL], color=INK, **dash)
    ax.plot([FR[0], FR[0] + 90], [cR, cR], color=INK, **dash)
    # walls: undeformed verticals and the inward shift of each top corner
    zt = h + b + 58
    ax.plot([FL[0], FL[0]], [cL, zt + 12], color=INK, **dash)
    ax.plot([FR[0], FR[0]], [cR, zt + 12], color=INK, **dash)
    ax.plot([TL[0], TL[0]], [h, zt + 12], color=INK, **dash)
    ax.plot([TR[0], TR[0]], [h + b, zt + 12], color=INK, **dash)
    dim((FL[0], zt), (TL[0], zt), "aL", dy=18)
    dim((FR[0], zt), (TR[0], zt), "aR", dy=18)
    # the sharp corners behind the rounding and the chamfer
    ax.plot([pr.A1[0], TL[0], pr.A2[0]], [pr.A1[1], TL[1], pr.A2[1]], color=INK, **dash)
    ax.plot([pr.B2[0], TR[0], pr.B1[0]], [pr.B2[1], TR[1], pr.B1[1]], color=INK, **dash)
    # height of the left top corner, top tilt, oxide
    ax.plot([TL[0], TR[0] + 78], [h, h], color=INK, **dash)
    xh = -60.0
    dim((xh, 0), (xh, h), "h", dx=14, ha="left", colour="white")
    dim((TR[0] + 62, h), (TR[0] + 62, h + b), "b", dx=14, ha="left")
    xt = 110.0
    dim((xt, h + b * (xt - TL[0]) / (TR[0] - TL[0])), (xt, H), "t", dx=14, ha="left")
    # rounding radius (left) and chamfer lengths (right)
    c = Profile._arc_center(TL, ((FL[0] - TL[0]) / pr.lenL, (FL[1] - TL[1]) / pr.lenL),
                            ((TR[0] - TL[0]) / pr.lenT, (TR[1] - TL[1]) / pr.lenT), pL, pr.thL)
    ang = math.radians(128)
    tip = (c[0] + pL * math.cos(ang), c[1] + pL * math.sin(ang))
    ax.annotate("", tip, c, arrowprops=dict(arrowstyle="->", color="white", lw=1.0, shrinkA=0, shrinkB=0))
    ax.scatter([c[0]], [c[1]], s=8, color="white", zorder=5)
    ax.text(c[0] + 12, c[1] - 4, "pL", ha="left", va="top", color="white", **mono)
    ax.text(0.5 * (TR[0] + pr.B2[0]), 0.5 * (TR[1] + pr.B2[1]) + 9, "qR", ha="center", va="bottom", color=INK, **mono)
    ax.set_xlim(X0, X1)
    ax.set_ylim(Z0, H + 50)
    ax.set_aspect("equal")
    ax.axis("off")
    save(fig, "mark")


def save(fig, name):
    fig.savefig(os.path.join(HERE, name + ".svg"), facecolor="white", metadata={"Date": None})
    plt.close(fig)


if __name__ == "__main__":
    figure_readings()
    finish_spectrum()
    figure_small()
    figure_mark()
