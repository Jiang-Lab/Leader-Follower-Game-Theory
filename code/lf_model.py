"""
Shared model code for the theory figures (Figs. 2-5, S1-S6).

Growth-coupled leader-follower game (main text Eqs. 1-5):
    f_L(x) = alpha*x + beta*(1-x)          leader fitness
    f_F(x) = gamma*x + delta*(1-x)         follower fitness
    dx/dt  = x(1-x) [f_L - f_F] (1 - N/K)
    dN/dt  = N [x f_L + (1-x) f_F] (1 - N/K)

Integrated with a vectorized fixed-step RK4 on (x, ln N), so that many
parameter sets (doses, grid points) are solved at once.

Treatment: a perturbation of payoff entry p with normalized dose d in [0, 1]
is p -> (1 - d) p (main text Eq. 9).
Objectives at time T (Eqs. 10-11):
    R_N = (N_untreated - N_treated) / N_untreated * 100   [%]
    S_L = 1 - x(T)
    J   = 0.5 * (R_N/100 + S_L)
"""
import numpy as np

# ---------------- default settings used throughout the theory figures ----------------
ALPHA, DELTA = 0.4, 0.6          # intrinsic growth rates (illustrative)
K = 1e10                          # carrying capacity
N0 = 1e5                          # initial population
X0 = 0.5                          # initial leader fraction

# representative (beta, gamma) for each regime
REGIMES = {
    "Snowdrift":          (0.8, 0.8),
    "Harmony":            (0.8, 0.2),
    "Coordination":       (0.2, 0.2),
    "Prisoner's Dilemma": (0.2, 0.8),
}
REGIME_ORDER = list(REGIMES)

# single-entry perturbations, in the order used in the figures
TARGETS = ["alpha", "beta", "gamma", "delta"]
TARGET_LABEL = {"alpha": r"$R_x\alpha$", "beta": r"$R_x\beta$",
                "gamma": r"$R_x\gamma$", "delta": r"$R_x\delta$"}
TARGET_COLOR = {"alpha": "#000000", "beta": "#0072B2", "gamma": "#009E73", "delta": "#D55E00"}
TARGET_STYLE = {"alpha": "-", "beta": "--", "gamma": "-.", "delta": ":"}


def simulate(alpha, beta, gamma, delta, x0=X0, n0=N0, k=K, T=15.0, dt=0.01,
             return_traj=False, n_save=None):
    """Integrate the growth-coupled model. All parameter arguments broadcast.

    Returns x(T), N(T) (arrays with the broadcast shape). With return_traj=True
    also returns t and the trajectories x(t), N(t) (time along axis 0).
    """
    a, b, g, d, x = np.broadcast_arrays(*[np.asarray(v, float) for v in
                                          (alpha, beta, gamma, delta, x0)])
    a, b, g, d = a.copy(), b.copy(), g.copy(), d.copy()
    x = x.copy()
    u = np.full_like(x, np.log(n0))          # u = ln N
    lnk = np.log(k)

    def rhs(x, u):
        h = 1.0 - np.exp(u - lnk)
        fl = a * x + b * (1.0 - x)
        ff = g * x + d * (1.0 - x)
        return x * (1.0 - x) * (fl - ff) * h, (x * fl + (1.0 - x) * ff) * h

    n_steps = int(round(T / dt))
    if return_traj:
        every = max(1, n_steps // n_save) if n_save else 1
        ts, xs, us = [0.0], [x.copy()], [u.copy()]
    for i in range(n_steps):
        k1x, k1u = rhs(x, u)
        k2x, k2u = rhs(x + 0.5 * dt * k1x, u + 0.5 * dt * k1u)
        k3x, k3u = rhs(x + 0.5 * dt * k2x, u + 0.5 * dt * k2u)
        k4x, k4u = rhs(x + dt * k3x, u + dt * k3u)
        x = np.clip(x + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x), 0.0, 1.0)
        u = np.minimum(u + dt / 6.0 * (k1u + 2 * k2u + 2 * k3u + k4u), lnk)
        if return_traj and ((i + 1) % every == 0 or i == n_steps - 1):
            ts.append((i + 1) * dt); xs.append(x.copy()); us.append(u.copy())
    if return_traj:
        return x, np.exp(u), np.array(ts), np.array(xs), np.exp(np.array(us))
    return x, np.exp(u)


def perturbed(beta, gamma, doses, alpha=ALPHA, delta=DELTA):
    """Payoff entries after applying doses = {'alpha': d1, 'beta': d2, ...} (broadcast)."""
    p = {"alpha": alpha, "beta": beta, "gamma": gamma, "delta": delta}
    out = {}
    for name, val in p.items():
        dose = doses.get(name, 0.0)
        out[name] = (1.0 - np.asarray(dose, float)) * val
    return out["alpha"], out["beta"], out["gamma"], out["delta"]


def treatment_outcome(beta, gamma, doses, T=15.0, alpha=ALPHA, delta=DELTA, x0=X0, dt=0.01):
    """R_N (%), S_L and J at time T for the given doses (broadcast)."""
    _, n_u = simulate(alpha, beta, gamma, delta, x0=x0, T=T, dt=dt)
    x_t, n_t = simulate(*perturbed(beta, gamma, doses, alpha, delta), x0=x0, T=T, dt=dt)
    r_n = (n_u - n_t) / n_u * 100.0
    s_l = 1.0 - x_t
    return r_n, s_l, 0.5 * (r_n / 100.0 + s_l), x_t


def regime_name(alpha, beta, gamma, delta):
    if beta > delta and gamma > alpha:
        return "Snowdrift"
    if beta > delta and gamma < alpha:
        return "Harmony"
    if beta < delta and gamma < alpha:
        return "Coordination"
    return "Prisoner's Dilemma"


def interior_equilibrium(alpha, beta, gamma, delta):
    den = (beta - delta) + (gamma - alpha)
    return np.nan if abs(den) < 1e-12 else (beta - delta) / den


def crossing_dose(target, alpha, beta, gamma, delta):
    """Regime-transition dose d* (main text Table 2); nan if no crossing in (0, 1)."""
    ratio = {"alpha": gamma / alpha, "gamma": alpha / gamma,
             "beta": delta / beta, "delta": beta / delta}[target]
    d = 1.0 - ratio
    return d if 0.0 < d < 1.0 else np.nan


def set_style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 12, "axes.labelsize": 15, "axes.titlesize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11,
        "axes.linewidth": 1.4, "lines.linewidth": 2.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def panel_letter(ax, letter, x=-0.18, y=1.08, size=18):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="left")


def save(fig, name, outdir=None):
    import os
    outdir = outdir or os.path.join("..", "figures")
    os.makedirs(outdir, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), bbox_inches="tight", dpi=300)
    print("Saved", os.path.join(outdir, name + ".png/.pdf"))


def composite_heatmap(fig, cell, r_n, s_l, best=None, xlabel="", ylabel="", doses=None,
                      show_ticks=True):
    """Two stacked dose-dose heatmaps in one grid cell: tumor reduction R_N (%) on top,
    leader suppression S_L (%) below. Arrays are indexed [dose_y, dose_x].
    best = (iy, ix) marks the dose pair maximizing J with an open red square in both maps."""
    import matplotlib.pyplot as plt
    doses = np.linspace(0, 1, r_n.shape[0]) if doses is None else doses
    step = doses[1] - doses[0]
    ext = [doses[0] - step / 2, doses[-1] + step / 2, doses[0] - step / 2, doses[-1] + step / 2]
    sub = cell.subgridspec(2, 1, hspace=0.0)
    axes = []
    for k, (arr, lab) in enumerate([(np.clip(r_n, 0, 100), "R_N"), (100 * s_l, "S_L")]):
        ax = fig.add_subplot(sub[k])
        im = ax.imshow(arr, origin="lower", extent=ext, cmap="magma", vmin=0, vmax=100,
                       interpolation="nearest", aspect="auto")
        if best is not None:
            iy, ix = best
            ax.plot(doses[ix], doses[iy], "s", ms=9, mfc="none", mec="red", mew=2.2)
        ax.set_xticks([0, 0.5, 1] if (show_ticks and k == 1) else [])
        ax.set_yticks([0, 1] if show_ticks else [])
        ax.tick_params(labelsize=8, length=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        axes.append(ax)
    axes[1].set_xlabel(xlabel, fontsize=11, labelpad=1)
    axes[0].set_ylabel(ylabel, fontsize=11, labelpad=1)
    axes[0].yaxis.set_label_coords(-0.22, 0.0)
    axes[0].axhline(ext[2], color="white", lw=1.5)
    return axes, im
