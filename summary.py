"""Every number quoted in README.md and pairs/README.md, recomputed from pairs/ and records/.

    python summary.py                 prints them, section by section
    python summary.py --catalogue     prints pairs/README.md
"""
import csv
import glob
import itertools
import json
import math
import os
import sys

import numpy as np

import model as PR

HERE = os.path.dirname(os.path.abspath(__file__))
CAPS = [10.0, 20.0, 40.0, 80.0]
GRID = 50          # nm, the wavelength step of records/band/
TWELVE = [500 + 45 * k for k in range(12)]   # a choice of twelve colours between 500 and 1000 nm
# what the sensor reads, as noise on (position, visibility, strength); None: not read, 0: exact
CASES = [("everything, exactly", (0.0, 0.0, 0.0)),
         ("positions to 0.01 nm, visibility and strength to 1e-4", (0.01, 1e-4, 1e-4)),
         ("positions to 0.01 nm, visibility and strength to 1e-3", (0.01, 1e-3, 1e-3)),
         ("positions to 0.01 nm, visibility and strength to 1e-2", (0.01, 1e-2, 1e-2)),
         ("positions only, exactly", (0.0, None, None)),
         ("positions only, to 0.01 nm", (0.01, None, None))]


def load():
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, "pairs", "*.txt"))):
        p = PR.read(f)
        p["name"] = os.path.basename(f)[:-4]
        p["record"] = json.load(open(os.path.join(HERE, "records", p["name"] + ".json")))
        out.append(p)
    return out


def covers(pair, cap, channels, case):
    """Does this pair show two indistinguishable marks for a sensor reading `channels` with noise
    `case`, when deformations up to `cap` are allowed?  (MODEL.md: fewer channels, more noise, an
    unread quantity and a larger cap all keep a pair indistinguishable.)"""
    if pair["cap"] > cap or not set(channels) <= set(pair["channels"]):
        return False
    for have, want in zip(pair["sigma"], case):
        if want is None:
            continue
        if have is None or have > want:
            return False
    return True


def best(pairs, cap, channels, case, floor=0.0):
    c = [p for p in pairs if covers(p, cap, channels, case) and p["separation"] > floor]
    return max(c, key=lambda p: p["separation"]) if c else None


def subsets():
    return [c for k in range(1, 9) for c in itertools.combinations(PR.NAMES, k)]


def band(name):
    rows = list(csv.DictReader(open(os.path.join(HERE, "records", "band", name + ".csv"))))
    for r in rows:
        for k in r:
            if k != "polarisation":
                r[k] = float(r[k])
    return rows


def band_difference(pair, rows, lo=500.0, hi=1000.0):
    """Differences between the two marks' readings over a band scan, in the channels the pair
    does not claim: (median, largest, where the largest is, number of channels)."""
    own = [(l, p) for l, p in PR.CHANNELS if "%d%s" % (l, p) in pair["channels"]]
    out = []
    for r in rows:
        if not lo <= r["wavelength_nm"] <= hi or (r["wavelength_nm"], r["polarisation"]) in own:
            continue
        s2 = -r["shift1_nm"] if pair["kind"] == "mirror" else r["shift2_nm"]
        d = abs(float(PR.wrap(pair["placement1"] + r["shift1_nm"] - pair["placement2"] - s2)))
        out.append((d, "%g %s" % (r["wavelength_nm"], r["polarisation"])))
    out.sort()
    return float(np.median([d for d, _ in out])), out[-1][0], out[-1][1], len(out)


def dots(pair):
    return "".join("●" if c in pair["channels"] else "○" for c in PR.NAMES)


def noise_words(sigma):
    if sigma == [0.0, 0.0, 0.0]:
        return "identical"
    if sigma == [0.0, None, None]:
        return "positions identical"
    if sigma == [0.01, 0.0, 0.0]:
        return "within 0.01 nm"
    return "within 0.01 nm, %g" % sigma[1]


def catalogue(pairs):
    out = []
    w = out.append
    general = [p["pair"] for p in pairs if p["kind"] != "mirror"]
    w("# The 51 pairs\n")
    w("One plain-text file per pair: two marks, their placements, the channels that cannot tell them apart and the noise "
      "that claim assumes. Definitions are in [MODEL.md](../MODEL.md); both solvers' data for each pair are in "
      "[records/](../records/) under the same name.\n")
    w("A **mirror** pair is a mark and its mirror image, which have the same visibility and strength in every channel. "
      "For identical data the mark must also shift every blind channel by the same amount: those equal-shift equations were "
      "solved by Newton's method on a rigorous solver, and the mirror image is then placed twice the shift away. "
      "For a mirror pair within noise the shifts need only agree to within the position noise: a nearly blind mark, scaled "
      "towards its symmetric form until the spread fits. Pairs %s are **general** pairs, two different marks: "
      "pair 08 solves for equal readings in all eight channels, pairs 06 and 07 are two nearby marks along the direction "
      "the data constrain least, and pair 39 minimises the difference in the data at a fixed distance.\n" % ", ".join(general))
    w("- **cap**: no deformation parameter of either mark exceeds it (nm).")
    w("- **blind channels**: filled where the pair claims the channel, in the order 532 TE, 532 TM, 633 TE, 633 TM, 775 TE, 775 TM, 850 TE, 850 TM.")
    w("- **data**: *identical*: every position, visibility and strength agrees. *positions identical*: the eight positions "
      "agree and the intensities are not compared. *within 0.01 nm*: `d <= 1` with that noise on every position and equal "
      "visibility and strength; a second number is the noise on visibility and strength.")
    w("- **apart**: distance between the two placements (nm). A lower bound on what those channels can hide: the pairs were built to pass 0.1 nm, not to be as far apart as possible.")
    w("- **centroids**: the same distance measured between the centroids of the two lines' cross-sections (nm).")
    w("- **imbalance**: the largest `|e+ - e-| / (e+ + e-)` over the blind channels. The two marks of a mirror pair have opposite imbalances.")
    w("- **other**: the largest difference between the two marks' readings in those of the eight channels the pair does not claim (nm).")
    w("- **median**, **largest**: the median and the largest difference between the two marks' readings at other wavelengths, "
      "every %d nm from 500 to 1000 nm in both polarisations (nm). The largest is usually where that line's first orders are weak.\n" % GRID)
    w("| pair | cap | blind channels | data | apart | centroids | imbalance | other | median | largest |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for p in pairs:
        r = p["record"]
        idx = [PR.NAMES.index(c) for c in p["channels"]]
        imb = float(np.abs(np.array(r["bie"]["imbalance"])[:, idx]).max())
        other = r["largest_reading_difference_outside_claim"]
        bd = band_difference(p, band(p["name"]))
        w("| [%s](%s.txt) | %g | %s | %s | %.4f | %.4f | %.1e | %s | %.3f | %.2f |" % (
            p["pair"], p["name"], p["cap"], dots(p), noise_words(p["sigma"]), p["separation"],
            r["separation_between_centroids"], imb, "" if other is None else "%.3f" % other, bd[0], bd[1]))
    return "\n".join(out) + "\n"


def main():
    pairs = load()
    if "--catalogue" in sys.argv:
        sys.stdout.write(catalogue(pairs))
        return
    P = {p["pair"]: p for p in pairs}
    head = P["01"]
    r = head["record"]
    m = head["mark1"]
    sh = np.array(r["bie"]["shift"])
    print("== the pair of the first figure (pair 01)")
    print("parameters of mark 1:", " ".join("%.10g" % v for v in m))
    floor = 0.5 * (m[3] + m[4])
    print("line height above the floor %.1f nm, width at the feet %.1f nm, oxide above the line %.1f nm" % (m[1] - floor, m[0], m[2]))
    print("wall angle from the wafer plane %.1f deg (walls lean %.1f and %.1f nm), shoulder radii %.1f and %.1f nm" % (
        math.degrees(math.atan((m[1] - floor) / (0.5 * (m[5] + m[6])))), m[5], m[6], m[8], m[9]))
    print("left minus right: wall %.4f, shoulder %.4f, floor %.4f, top tilt %.4f nm" % (m[5] - m[6], m[8] - m[9], m[3] - m[4], m[7]))
    print("placements %.10f and %.10f, apart %.10f nm" % (head["placement1"], head["placement2"], head["separation"]))
    print("shift of mark 1 per channel (nm):", " ".join("%.10f" % v for v in sh[0]))
    print("mean shift %.10f nm; spread over the eight channels %.1e nm (mark 1), %.1e nm (mark 2)" % (
        sh[0].mean(), sh[0].max() - sh[0].min(), sh[1].max() - sh[1].min()))
    X = sh + np.array([[head["placement1"]], [head["placement2"]]])
    print("largest difference between the two marks: reading %.1e nm, visibility %.1e, strength %.1e" % (
        np.abs(X[0] - X[1]).max(), np.abs(np.array(r["bie"]["visibility"][0]) - np.array(r["bie"]["visibility"][1])).max(),
        np.abs(np.array(r["bie"]["strength"][0]) - np.array(r["bie"]["strength"][1])).max()))
    print("visibility:", " ".join("%.10f" % v for v in r["bie"]["visibility"][0]))
    print("strength:  ", " ".join("%.4f" % v for v in r["bie"]["strength"][0]))
    print("imbalance mark 1:", " ".join("%+.2e" % v for v in r["bie"]["imbalance"][0]))
    print("imbalance mark 2:", " ".join("%+.2e" % v for v in r["bie"]["imbalance"][1]))
    imb = np.abs(np.array(r["bie"]["imbalance"][0]))
    print("imbalance magnitude from %.1e (%s) to %.1e (%s)" % (imb.min(), PR.NAMES[int(imb.argmin())], imb.max(), PR.NAMES[int(imb.argmax())]))
    print("separation between centroids %.4f nm (centroid offsets %+.4f, %+.4f)" % (
        r["separation_between_centroids"], r["centroid_offset"][0], r["centroid_offset"][1]))
    print("two solvers on pair 01: readings within %.1e nm (order 10)" % r["agreement_bie_fem10"]["shift"])

    print("\n== one asymmetry at a time (records/asymmetries.json)")
    for c in json.load(open(os.path.join(HERE, "records", "asymmetries.json")))["cases"]:
        s = np.array(c["shift"])
        print("%-10s from %+.4f to %+.4f nm, spread %.4f   (%s)" % (c["name"], s.min(), s.max(), s.max() - s.min(), c["what"]))

    p19 = P["19"]
    q = p19["mark1"]
    r19 = p19["record"]
    X19 = np.array(r19["bie"]["shift"]) + np.array([[p19["placement1"]], [p19["placement2"]]])
    out = [j for j, c in enumerate(PR.NAMES) if c not in p19["channels"]]
    print("\n== pair 19, the nearly perfect line")
    print("walls lean %.2f and %.2f nm, shoulders rounded by %.2f and %.2f nm, floor tilt %.3f and %.3f, top tilt %.3f; largest deformation %.2f nm" % (
        q[5], q[6], q[8], q[9], q[3], q[4], q[7], max(np.abs(p19["mark1"][3:]).max(), np.abs(p19["mark2"][3:]).max())))
    print("blind channels: %s; apart %.4f nm; common reading %.4f nm" % (" ".join(p19["channels"]), p19["separation"],
                                                                    X19[0, [j for j in range(8) if j not in out]].mean()))
    print("the other channels differ by: %s" % ", ".join("%s %.4f nm" % (PR.NAMES[j], abs(X19[0, j] - X19[1, j])) for j in out))
    p08 = P["08"]
    print("pair 08 (positions only): largest deformation %.2f nm, apart %.4f nm" % (
        max(np.abs(p08["mark1"][3:]).max(), np.abs(p08["mark2"][3:]).max()), p08["separation"]))

    print("\n== all 51 pairs")
    print("pairs %d: mirror %d, general %d" % (len(pairs), sum(p["kind"] == "mirror" for p in pairs), sum(p["kind"] != "mirror" for p in pairs)))
    print("identical data %d, positions identical %d, within noise %d" % (
        sum(p["sigma"] == [0.0, 0.0, 0.0] for p in pairs), sum(p["sigma"] == [0.0, None, None] for p in pairs),
        sum(p["sigma"][0] == 0.01 for p in pairs)))
    print("more than 0.1 nm apart: %d; on all eight channels: %d, of which more than 0.1 nm apart: %d" % (
        sum(p["separation"] > 0.1 for p in pairs), sum(len(p["channels"]) == 8 for p in pairs),
        sum(len(p["channels"]) == 8 and p["separation"] > 0.1 for p in pairs)))
    print("largest separation %.4f nm (pair %s)" % max((p["separation"], p["pair"]) for p in pairs))
    print("two solvers, all pairs: readings within %.1e nm, visibility %.1e, strength %.1e (finite elements order 10, claimed channels)" % tuple(
        max(p["record"]["agreement_bie_fem10"][k] for p in pairs) for k in ("shift", "visibility", "strength")))
    if all("agreement_bie_fem8" in p["record"] for p in pairs):
        print("                        readings within %.1e nm, visibility %.1e, strength %.1e (finite elements order 8, all eight channels)" % tuple(
            max(p["record"]["agreement_bie_fem8"][k] for p in pairs) for k in ("shift", "visibility", "strength")))
    print("boundary integrals, two resolutions: readings within %.1e nm, visibility %.1e, strength %.1e" % tuple(
        max(p["record"]["bie"]["change_with_resolution"][k] for p in pairs) for k in ("shift", "visibility", "strength")))
    eq = [max(p["record"]["largest_difference_where_equal"].values()) for p in pairs if p["record"]["largest_difference_where_equal"]]
    print("quantities claimed equal differ by at most %.1e (%d pairs)" % (max(eq), len(eq)))
    d2 = [p["record"]["d2"] for p in pairs if p["record"]["d2"] is not None]
    print("d^2 of the within-noise pairs: %.4f to %.4f (%d pairs)" % (min(d2), max(d2), len(d2)))
    print("apart by more than 0.1 nm between centroids: %d of %d" % (sum(p["record"]["separation_between_centroids"] > 0.1 for p in pairs), len(pairs)))

    print("\n== all eight channels: largest separation of a pair, by what is read and by cap (nm)")
    for label, case in CASES:
        row = []
        for cap in CAPS:
            b = best(pairs, cap, PR.NAMES, case)
            row.append("%.3f (%s)" % (b["separation"], b["pair"]) if b else "none")
        print("%-58s %s" % (label, " | ".join(row)))

    print("\n== everything read exactly: channel sets with a pair more than 0.1 nm apart, by size and cap")
    subs = subsets()
    for k in range(1, 9):
        row = []
        for cap in CAPS:
            ss = [c for c in subs if len(c) == k]
            row.append("%d of %d" % (sum(best(pairs, cap, c, CASES[0][1], 0.1) is not None for c in ss), len(ss)))
        print("%d channels: %s" % (k, " | ".join(row)))
    te = ("532TE", "633TE", "775TE", "850TE")
    tm = ("532TM", "633TM", "775TM", "850TM")
    for lab, c in (("four colours, TE only", te), ("four colours, TM only", tm)):
        b = best(pairs, 10.0, c, CASES[0][1], 0.1)
        print("%s, cap 10: %s" % (lab, "%.3f nm (pair %s)" % (b["separation"], b["pair"]) if b else "none"))
    tot = dec = 0
    open_ = {}
    for label, case in CASES:
        for cap in CAPS:
            for c in subs:
                tot += 1
                if best(pairs, cap, c, case, 0.1) is not None:
                    dec += 1
                else:
                    open_.setdefault((label, cap), []).append(c)
    print("\n== every cap, channel set and reading: %d of %d have a pair more than 0.1 nm apart; %d do not" % (dec, tot, tot - dec))
    for (label, cap), cs in open_.items():
        print("   %s, cap %g: %d sets without one (sizes %s)" % (label, cap, len(cs), sorted(set(len(c) for c in cs))))

    rows = list(csv.DictReader(open(os.path.join(HERE, "records", "spectrum-01.csv"))))
    lam = np.array([float(x["wavelength_nm"]) for x in rows])
    pol = np.array([x["polarisation"] for x in rows])
    s = np.array([float(x["shift1_nm"]) for x in rows])
    eff = np.array([float(x["eplus1"]) + float(x["eminus1"]) for x in rows])
    d0 = sh[0].mean()

    print("\n== more colours (records/spectrum-01.csv: mark 1 of pair 01, 1 nm steps; the reading is measured from the mark)")
    band_ = (lam >= 500) & (lam <= 1000)
    print("every wavelength from 500 to 1000 nm, both polarisations: %d readings from %.4f to %.4f nm, mean %.4f; no two differ by more than %.4f nm" % (
        band_.sum(), s[band_].min(), s[band_].max(), s[band_].mean(), s[band_].max() - s[band_].min()))
    print("weakest first orders in that range: e+ + e- = %.3f" % eff[band_].min())
    twelve = TWELVE
    k12 = np.isin(lam, twelve)
    x12 = s[k12]
    print("a choice of twelve colours, %d to %d nm in steps of %d nm, both polarisations: %d readings from %.4f to %.4f nm," % (
        twelve[0], twelve[-1], twelve[1] - twelve[0], x12.size, x12.min(), x12.max()))
    print("   spread %.4f nm, mean %.4f nm, standard deviation %.4f nm; the mirror image differs from the mark in one of them by at most %.4f nm" % (
        x12.max() - x12.min(), x12.mean(), x12.std(), 2 * np.abs(x12 - d0).max()))
    print("   anywhere from 500 to 1000 nm the mirror image differs from the mark by at most %.4f nm" % (2 * np.abs(s[band_] - d0).max()))
    fem = list(csv.DictReader(open(os.path.join(HERE, "records", "spectrum-01-fem.csv"))))
    bie_ = {(float(x["wavelength_nm"]), x["polarisation"]): float(x["shift1_nm"]) for x in rows}
    dif = [abs(float(x["shift1_nm"]) - bie_[(float(x["wavelength_nm"]), x["polarisation"])]) for x in fem]
    lamf = sorted(set(float(x["wavelength_nm"]) for x in fem))
    print("two solvers on this spectrum: finite elements (order 8) every %g nm from %g to %g nm, %d readings, within %.1e nm of the boundary integrals" % (
        lamf[1] - lamf[0], lamf[0], lamf[-1], len(fem), max(dif)))

    print("\n== outside that range (records/spectrum-01.csv)")
    for lo, hi in ((400, 1100), (500, 1100), (520, 1100), (500, 1000), (520, 1000)):
        k = (lam >= lo) & (lam <= hi)
        i0, i1 = np.argmin(np.where(k, s, np.inf)), np.argmax(np.where(k, s, -np.inf))
        print("%d to %d nm: shift from %+.4f (%g %s) to %+.4f (%g %s); farthest from %.4f by %.4f; two marks differ by up to %.4f nm" % (
            lo, hi, s[i0], lam[i0], pol[i0], s[i1], lam[i1], pol[i1], d0, np.abs(s[k] - d0).max(), 2 * np.abs(s[k] - d0).max()))
    i = np.argmin(eff)
    k633 = [float(x["eplus1"]) + float(x["eminus1"]) for x in rows if float(x["wavelength_nm"]) == 633.0]
    print("weakest first orders: e+ + e- = %.2e at %g nm %s (at 633 nm: %.3f TE, %.3f TM)" % (eff[i], lam[i], pol[i], k633[0], k633[1]))
    from bie.materials import n_sio2
    lam0 = lam[i]
    print("half a wavelength deep: 2 x n_oxide(%g nm) x line height = 2 x %.4f x %.2f = %.1f nm" % (
        lam0, n_sio2(lam0), m[1] - floor, 2 * n_sio2(lam0) * (m[1] - floor)))
    k532 = [float(x["eplus1"]) + float(x["eminus1"]) for x in rows if float(x["wavelength_nm"]) == 532.0]
    print("first orders at the weakest point against 532 nm: 1/%.0f (TE), 1/%.0f (TM)" % (k532[0] / eff[i], k532[1] / eff[i]))

    print("\n== all pairs at other wavelengths (records/band/, every %d nm from 500 to 1000 nm, channels outside each pair's claim)" % GRID)
    bd = [(band_difference(p, band(p["name"])), p["pair"]) for p in pairs]
    med = sorted((b[0], n) for b, n in bd)
    print("median difference between the two marks' readings: from %.3f nm (pair %s) to %.3f nm (pair %s); middle pair %.3f nm" % (
        med[0][0], med[0][1], med[-1][0], med[-1][1], med[len(med) // 2][0]))
    print("pairs whose median difference exceeds 0.01 nm: %d of %d; 0.05 nm: %d" % (
        sum(m > 0.01 for m, _ in med), len(med), sum(m > 0.05 for m, _ in med)))
    big = sorted((b[1], n) for b, n in bd)
    print("largest difference: from %.3f nm (pair %s) to %.2f nm (pair %s)" % (big[0][0], big[0][1], big[-1][0], big[-1][1]))
    print("channels compared per pair: %d to %d" % (min(b[3] for b, _ in bd), max(b[3] for b, _ in bd)))

    print("\n== order imbalance, mirror pairs (largest over the blind channels)")
    im = []
    for p in pairs:
        if p["kind"] != "mirror":
            continue
        idx = [PR.NAMES.index(c) for c in p["channels"]]
        im.append((float(np.abs(np.array(p["record"]["bie"]["imbalance"])[0, idx]).max()), p["pair"]))
    print("from %.1e (pair %s) to %.1e (pair %s); below 1e-4 in %d of %d" % (min(im)[0], min(im)[1], max(im)[0], max(im)[1],
                                                                          sum(v < 1e-4 for v, _ in im), len(im)))


if __name__ == "__main__":
    main()
