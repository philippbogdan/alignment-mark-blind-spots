"""Panel quadrature tables: Gauss-Legendre nodes, log-product weights, interpolation matrices."""
import math
import numpy as np


def gl(n):
    t, w = np.polynomial.legendre.leggauss(n)
    return t, w


def legendre_vander(t, n):
    return np.polynomial.legendre.legvander(t, n - 1)    # (len(t), n)


def log_moments(a, n):
    """m_k = int_{-1}^{1} P_k(t) log|t - a| dt, k < n, for a in (-1, 1)."""
    # Legendre functions of the second kind on the cut, forward recurrence
    Q = np.zeros(n + 1)
    Q[0] = 0.5 * math.log((1 + a) / (1 - a))
    if n >= 1:
        Q[1] = a * Q[0] - 1.0
    for l in range(1, n):
        Q[l + 1] = ((2 * l + 1) * a * Q[l] - l * Q[l - 1]) / (l + 1)
    m = np.zeros(n)
    m[0] = (1 - a) * math.log(1 - a) + (1 + a) * math.log(1 + a) - 2.0
    for k in range(1, n):
        m[k] = 2.0 / (2 * k + 1) * (Q[k + 1] - Q[k - 1])
    return m


def log_weights(t):
    """W[i, j]: sum_j W[i,j] f(t_j) = int_{-1}^1 f(t) log|t - t_i| dt for polynomials f of degree < n."""
    n = len(t)
    V = legendre_vander(t, n).T          # V[k, j] = P_k(t_j)
    W = np.empty((n, n))
    for i in range(n):
        W[i] = np.linalg.solve(V, log_moments(t[i], n))
    return W


def interp_matrix(t_from, t_to):
    """Lagrange interpolation (barycentric) from nodes t_from to points t_to."""
    t_from = np.asarray(t_from)
    t_to = np.asarray(t_to)
    n = len(t_from)
    wb = np.array([1.0 / np.prod(t_from[j] - np.delete(t_from, j)) for j in range(n)])
    D = t_to[:, None] - t_from[None, :]
    exact = D == 0
    D = np.where(exact, 1.0, D)
    C = wb[None, :] / D
    C = C / C.sum(axis=1, keepdims=True)
    rows = exact.any(axis=1)
    C[rows] = exact[rows].astype(float)
    return C


class Tables:
    def __init__(self, n=16, n_aux=16, aux_levels=9):
        self.n = n
        self.t, self.w = gl(n)
        self.WL = log_weights(self.t)
        self.logdt = np.log(np.abs(self.t[:, None] - self.t[None, :]) + np.eye(n))   # diag 0
        # children of a panel (halves): t in [-1,0] and [0,1]
        tc = np.concatenate([(self.t - 1) / 2, (self.t + 1) / 2])
        self.Pchild = interp_matrix(self.t, tc)          # (2n, n)
        # auxiliary graded rule toward t = -1 (u = (t+1)/2 in [0,1])
        ta, wa = gl(n_aux)
        edges = [0.0] + [2.0 ** (-k) for k in range(aux_levels, 0, -1)] + [1.0]
        us, ws = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            us.append(a + (b - a) * (ta + 1) / 2)
            ws.append((b - a) / 2 * wa)
        u = np.concatenate(us)
        wu = np.concatenate(ws)
        self.aux_t_minus = 2 * u - 1               # graded toward -1
        self.aux_w = 2 * wu                        # weights in t (sum = 2)
        self.aux_t_plus = -self.aux_t_minus[::-1]  # graded toward +1
        self.aux_w_plus = self.aux_w[::-1].copy()
        self.I_minus = interp_matrix(self.t, self.aux_t_minus)
        self.I_plus = interp_matrix(self.t, self.aux_t_plus)
