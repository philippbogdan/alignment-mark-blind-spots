"""Optical constants (own readers of the refractiveindex.info files in data/).

Si: tabulated n, k (Green 2008), linearly interpolated in wavelength, eps = (n + i k)^2.
SiO2: Sellmeier formula 1 (Malitson 1965).  Air: n = 1.  Time dependence exp(-i omega t).
"""
import os
import numpy as np

_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


def _read_block(path, key):
    """Return the lines following `key` inside the DATA block of a refractiveindex.info yml."""
    lines = open(path).read().splitlines()
    out, on = [], False
    for ln in lines:
        s = ln.strip()
        if s.startswith(key):
            on = True
            rest = s[len(key):].strip()
            if rest and rest != "|":
                out.append(rest)
            continue
        if on:
            if not s or (":" in s and not s[0].isdigit()):
                break
            out.append(s)
    return out


def _si_table():
    rows = [r.split() for r in _read_block(os.path.join(_DATA, "si-green-2008.yml"), "data:")]
    arr = np.array([[float(v) for v in r] for r in rows if len(r) == 3])
    return arr[:, 0] * 1000.0, arr[:, 1], arr[:, 2]   # wavelength in nm, n, k


_SI = _si_table()


def _sio2_coeffs():
    lines = open(os.path.join(_DATA, "sio2-malitson-1965.yml")).read().splitlines()
    for ln in lines:
        s = ln.strip()
        if s.startswith("coefficients:"):
            return np.array([float(v) for v in s.split(":", 1)[1].split()])
    raise RuntimeError("no Sellmeier coefficients")


_SIO2 = _sio2_coeffs()


def n_si(lam_nm):
    lam, n, k = _SI
    lam_nm = np.asarray(lam_nm, float)
    if np.any(lam_nm < lam[0]) or np.any(lam_nm > lam[-1]):
        raise ValueError("wavelength outside Si table")
    return np.interp(lam_nm, lam, n) + 1j * np.interp(lam_nm, lam, k)


def eps_si(lam_nm):
    return n_si(lam_nm) ** 2


def n_sio2(lam_nm):
    c = _SIO2
    l2 = (np.asarray(lam_nm, float) / 1000.0) ** 2
    n2 = 1.0 + c[0]
    for i in range(1, len(c) - 1, 2):
        n2 = n2 + c[i] * l2 / (l2 - c[i + 1] ** 2)
    return np.sqrt(n2)


def eps_sio2(lam_nm):
    return n_sio2(lam_nm) ** 2 + 0j
