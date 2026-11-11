"""Optical constants exactly as fixed in MODEL.md.

Silicon: Green 2008 table, n and k linearly interpolated in wavelength.
Silicon dioxide: Malitson 1965 Sellmeier formula (refractiveindex.info "formula 1").
Air: n = 1.  Permittivity eps = (n + i k)^2, time dependence exp(-i omega t).
"""
import os
import numpy as np

_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


def _read_si_table(path=os.path.join(_DATA, "si-green-2008.yml")):
    rows = []
    in_data = False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("data: |"):
                in_data = True
                continue
            if in_data:
                parts = s.split()
                if len(parts) == 3:
                    try:
                        rows.append([float(p) for p in parts])
                        continue
                    except ValueError:
                        pass
                if rows:
                    break
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2]   # wavelength (um), n, k


def _read_sellmeier(path=os.path.join(_DATA, "sio2-malitson-1965.yml")):
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("coefficients:"):
                return [float(v) for v in s.split(":", 1)[1].split()]
    raise RuntimeError("no Sellmeier coefficients")


_SI_WL, _SI_N, _SI_K = _read_si_table()
_SIO2_C = _read_sellmeier()


def n_si(lam_nm):
    """Complex refractive index of Si at wavelength lam_nm (nm), linear interpolation in wavelength."""
    lam_um = np.asarray(lam_nm, dtype=float) / 1000.0
    if np.any(lam_um < _SI_WL[0]) or np.any(lam_um > _SI_WL[-1]):
        raise ValueError("wavelength outside Si table")
    n = np.interp(lam_um, _SI_WL, _SI_N)
    k = np.interp(lam_um, _SI_WL, _SI_K)
    return n + 1j * k


def n_sio2(lam_nm):
    """Refractive index of SiO2 (Malitson Sellmeier, refractiveindex.info formula 1)."""
    lam_um = np.asarray(lam_nm, dtype=float) / 1000.0
    c = _SIO2_C
    n2 = 1.0 + c[0]
    l2 = lam_um ** 2
    for i in range(1, len(c), 2):
        n2 = n2 + c[i] * l2 / (l2 - c[i + 1] ** 2)
    return np.sqrt(n2)


def eps_si(lam_nm):
    return n_si(lam_nm) ** 2


def eps_sio2(lam_nm):
    return n_sio2(lam_nm) ** 2 + 0j


CHANNELS = [(532.0, "TE"), (532.0, "TM"), (633.0, "TE"), (633.0, "TM"),
            (775.0, "TE"), (775.0, "TM"), (850.0, "TE"), (850.0, "TM")]
CHANNEL_NAMES = ["%d%s" % (int(l), p) for (l, p) in CHANNELS]

if __name__ == "__main__":
    for lam in (532, 633, 775, 850):
        print(lam, n_si(lam), eps_si(lam), n_sio2(lam), eps_sio2(lam))
