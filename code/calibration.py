"""
Calibration of the leader-follower game model to the H1299 growth data
(Python version of the MATLAB scripts Fig5_calibration_monoculture_and_coculture_fits.m
and Table3_S2_parameters_both_calibrations.m; same data, model, bounds, weighting,
bootstrap and AICc).

Two calibrations:
    start = 0  including the first day: days 0-4, every series normalized to day 0 (main fit)
    start = 1  excluding the first day: days 1-4, every series normalized to day 1

Monocultures   y = y0 exp(r t), r in [0, 1], y0 within the data range, fitted to the
               replicate mean; replicate bootstrap (resample the 5 replicates).
Co-cultures    joint fit of all six co-cultures to the leader fraction x and N/N0, residuals
               divided by the pooled SD of each variable; alpha and delta fixed from the
               monocultures; beta, gamma in [-2, 3]; four cases
               H0 (beta, gamma shared), Hbeta (only beta differs), Hgamma (only gamma
               differs), Hbg (both differ); residual bootstrap; AICc.

The bootstrap uses NumPy's random generator, so bootstrap means and CIs agree with the
MATLAB output up to resampling noise (second or third decimal); best-fit values agree.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

# ----------------------------------------------------------------------------- data
# Monocultures normalized to day 0; rows = days 0-4, columns = 5 replicates
MONO = {
    ("Leaders", 1):   [[1, 1, 1, 1, 1],
                       [2.3510, 2.4155, 1.1853, 1.2026, 2.1348],
                       [3.0204, 3.2785, 2.8707, 3.0000, 2.8090],
                       [4.5306, 5.9361, 3.6293, 4.2716, 3.5955],
                       [6.8163, 6.4840, 5.5603, 5.9052, 4.4944]],
    ("Leaders", 2):   [[1, 1, 1, 1, 1],
                       [0.8143, 0.7173, 1.7026, 1.7974, 1.4607],
                       [1.4051, 0.8342, 2.2586, 2.3147, 1.6854],
                       [0.6257, 0.9541, 2.0560, 2.3534, 1.7978],
                       [1.0212, 1.0425, 2.7026, 2.4828, 1.9101]],
    ("Followers", 1): [[1, 1, 1, 1, 1],
                       [1.2577, 1.5942, 1.5055, 1.5240, 1.8333],
                       [1.7298, 1.0526, 1.2362, 1.1107, 2.7500],
                       [3.5769, 2.8321, 3.8376, 3.8376, 4.1667],
                       [6.6000, 4.9861, 6.7159, 6.3100, 6.6667]],
    ("Followers", 2): [[1, 1, 1, 1, 1],
                       [1.7077, 0.9530, 0.9041, 0.9041, 1.7500],
                       [1.3154, 0.9372, 1.2362, 1.2694, 2.6667],
                       [4.0192, 3.0511, 3.2657, 3.2030, 3.8333],
                       [4.7731, 3.5365, 6.2362, 6.4945, 6.5000]],
}
MONO = {k: np.array(v, float) for k, v in MONO.items()}

# Co-cultures: condition (1 untreated, 2 +TGF-beta1), ratio (1 = 1:1, 2 = 1:2, 3 = 1:5 L:F),
# leader counts, total counts on days 0-4 (one replicate per ratio)
COCULTURE = [
    (1, 1, [4.5e5, 7.6e5, 9.0e5, 1.0e6, 9.4e5], [8.7e5, 1.7e6, 2.7e6, 3.4e6, 4.3e6]),
    (1, 2, [3.6e5, 6.2e5, 5.8e5, 6.5e5, 6.4e5], [1.0e6, 2.0e6, 2.8e6, 3.8e6, 5.3e6]),
    (1, 3, [1.8e5, 2.7e5, 2.6e5, 3.3e5, 4.0e5], [9.9e5, 1.9e6, 3.0e6, 4.4e6, 7.3e6]),
    (2, 1, [4.5e5, 8.1e5, 7.8e5, 6.6e5, 5.3e5], [8.7e5, 1.9e6, 2.5e6, 2.8e6, 3.5e6]),
    (2, 2, [3.6e5, 5.4e5, 5.5e5, 5.5e5, 4.3e5], [1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.7e6]),
    (2, 3, [1.8e5, 2.6e5, 2.4e5, 2.7e5, 2.3e5], [9.9e5, 2.0e6, 3.0e6, 4.4e6, 6.7e6]),
]
CASES = ["H0", "Hbeta", "Hgamma", "Hbg"]
CASE_K = {"H0": 2, "Hbeta": 3, "Hgamma": 3, "Hbg": 4}
CASE_TEX = {"H0": r"$H_0$", "Hbeta": r"$H_\beta$", "Hgamma": r"$H_\gamma$", "Hbg": r"$H_{\beta\gamma}$"}
COND = {1: "untreated", 2: "+TGF-beta1"}


# ----------------------------------------------------------------------------- monocultures
def mono_data(start):
    """Monoculture replicates from day `start`, each normalized to its value on that day."""
    out = {}
    for k, T in MONO.items():
        T = T[start:]
        out[k] = T / T[0]
    return out


def mono_fit(T, n_boot, rng, t_grid=None):
    """Replicate bootstrap of y = y0 exp(r t) fitted to the replicate mean.
    Returns dict(r_mean, r_ci, fit, band) (fit/band on t_grid if given)."""
    t = np.arange(T.shape[0], dtype=float)
    lo, hi = T.min(), T.max()
    n_rep = T.shape[1]
    P = []
    for _ in range(n_boot):
        y = T[:, rng.integers(0, n_rep, n_rep)].mean(1)
        s = least_squares(lambda p: p[1] * np.exp(p[0] * t) - y, [0.5, np.clip(T[0, 0], lo, hi)],
                          bounds=([0, lo], [1, hi]), xtol=1e-12, ftol=1e-12, gtol=1e-12)
        P.append(s.x)
    P = np.array(P)
    out = dict(r_mean=P[:, 0].mean(), r_ci=np.percentile(P[:, 0], [2.5, 97.5]), P=P)
    if t_grid is not None:
        Y = P[:, 1:2] * np.exp(P[:, 0:1] * t_grid[None, :])
        out["fit"] = Y.mean(0)
        out["band"] = np.percentile(Y, [2.5, 97.5], 0)
    return out


def fit_monocultures(start, n_boot, rng, t_grid=None):
    """alpha (untreated, +TGF), delta (shared = untreated estimate) and the full fits."""
    M = {k: mono_fit(T, n_boot, rng, t_grid) for k, T in mono_data(start).items()}
    alpha = np.array([M[("Leaders", 1)]["r_mean"], M[("Leaders", 2)]["r_mean"]])
    delta = np.array([M[("Followers", 1)]["r_mean"]] * 2)
    return alpha, delta, M


# ----------------------------------------------------------------------------- co-cultures
def coculture_data(start, alpha, delta):
    """Series from day `start`: x = L/N and n = N/N(start); pooled SDs W for the weighting."""
    S = []
    for cond, ratio, L, N in COCULTURE:
        L = np.array(L, float)[start:]
        N = np.array(N, float)[start:]
        S.append(dict(cond=cond, ratio=ratio, x=L / N, n=N / N[0],
                      alpha=alpha[cond - 1], delta=delta[cond - 1]))
    allx = np.concatenate([s["x"][1:] for s in S])
    alln = np.concatenate([s["n"][1:] for s in S])
    return S, np.array([allx.std(), alln.std()])


def simulate(beta, gamma, alpha, delta, x0, t):
    """Leader fraction x(t) and N(t)/N(0) of the game model (no carrying-capacity term:
    N/K < 1e-3 over the assay)."""
    def rhs(_, z):
        x, n = z
        fl = alpha * x + beta * (1 - x)
        ff = gamma * x + delta * (1 - x)
        return [(fl - ff) * x * (1 - x), (x * fl + (1 - x) * ff) * n]
    s = solve_ivp(rhs, (t[0], t[-1]), [x0, 1.0], t_eval=t, rtol=1e-9, atol=1e-11)
    return s.y[0], s.y[1]


def fitness_difference(beta, gamma, alpha, delta, x):
    """f_L - f_F = (beta - delta) + (alpha - beta - gamma + delta) x."""
    return (alpha * x + beta * (1 - x)) - (gamma * x + delta * (1 - x))


def expand(p, case):
    """[beta_U, gamma_U, beta_T, gamma_T] for each case."""
    if case == "H0":
        return np.array([p[0], p[1], p[0], p[1]])
    if case == "Hbeta":
        return np.array([p[0], p[1], p[2], p[1]])
    if case == "Hgamma":
        return np.array([p[0], p[1], p[0], p[2]])
    return np.array(p[:4])


def residuals(q, S, W):
    t = np.arange(len(S[0]["x"]), dtype=float)
    r = []
    for s in S:
        b, g = q[2 * s["cond"] - 2], q[2 * s["cond"] - 1]
        x, n = simulate(b, g, s["alpha"], s["delta"], s["x"][0], t)
        r += [(x[1:] - s["x"][1:]) / W[0], (n[1:] - s["n"][1:]) / W[1]]
    return np.concatenate(r)


def fit_case(case, S, W, p0=None):
    """Best fit of one case; returns q = [beta_U, gamma_U, beta_T, gamma_T] and SSE."""
    k = CASE_K[case]
    p0 = 0.3 * np.ones(k) if p0 is None else p0
    s = least_squares(lambda p: residuals(expand(p, case), S, W), p0,
                      bounds=([-2] * k, [3] * k), xtol=1e-10, ftol=1e-10, gtol=1e-10)
    return expand(s.x, case), 2 * s.cost


def aicc(sse, n, k):
    kk = k + 1                         # fitted parameters + error variance
    return n * np.log(sse / n) + 2 * kk + 2 * kk * (kk + 1) / (n - kk - 1)


def boot_resample(S, q, c, rng):
    """Residual bootstrap: centered, inflated residuals resampled and added to the fit."""
    t = np.arange(len(S[0]["x"]), dtype=float)
    m = len(t) - 1
    Sb = []
    for s in S:
        b, g = q[2 * s["cond"] - 2], q[2 * s["cond"] - 1]
        x, n = simulate(b, g, s["alpha"], s["delta"], s["x"][0], t)
        rx = s["x"][1:] - x[1:]; rx = (rx - rx.mean()) * c
        rn = s["n"][1:] - n[1:]; rn = (rn - rn.mean()) * c
        Sb.append(dict(s, x=np.r_[s["x"][0], x[1:] + rx[rng.integers(0, m, m)]],
                       n=np.r_[1.0, n[1:] + rn[rng.integers(0, m, m)]]))
    return Sb


def bootstrap_case(case, q, S, W, n_boot, seed):
    """Residual-bootstrap parameter draws (n_boot x 4) for one case."""
    rng = np.random.default_rng(seed)
    n = 2 * (len(S[0]["x"]) - 1) * len(S)
    c = np.sqrt(n / (n - CASE_K[case]))
    P = []
    for _ in range(n_boot):
        try:
            P.append(fit_case(case, boot_resample(S, q, c, rng), W)[0])
        except Exception:
            pass
    return np.array(P)


# ----------------------------------------------------------------------------- game analysis
def xstar(alpha, beta, gamma, delta):
    """Interior equilibrium x* = (beta - delta) / ((beta - delta) + (gamma - alpha))."""
    den = (beta - delta) + (gamma - alpha)
    return np.nan if abs(den) < 1e-12 else (beta - delta) / den


def game(alpha, beta, gamma, delta):
    if beta > delta and gamma > alpha:
        return "Snowdrift"
    if beta > delta and gamma < alpha:
        return "Harmony"
    if beta < delta and gamma < alpha:
        return "Coordination"
    return "Prisoner's Dilemma"


def outcome(g, xs, x0max):
    """Stable state reached from the observed initial leader fractions."""
    if g == "Prisoner's Dilemma":
        return "x -> 0 (followers only, from any x0)"
    if g == "Harmony":
        return "x -> 1 (leaders only, from any x0)"
    if g == "Snowdrift":
        return f"x -> x* = {xs:.2f} (coexistence)"
    return (f"x -> 0 (followers only: all x0 <= {x0max:.2f} < x* = {xs:.2f})" if x0max < xs
            else f"depends on x0 (x* = {xs:.2f})")


def bootstrap_job(args):
    """Worker for parallel bootstraps: args = (case, q, S, W, n_boot, seed).
    Defined here (not in the notebook) so that it also works with the 'spawn' start method
    used on Windows and macOS."""
    case, q, S, W, n_boot, seed = args
    return bootstrap_case(case, q, S, W, n_boot, seed)
