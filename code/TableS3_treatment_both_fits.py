import numpy as np
import pandas as pd
from itertools import combinations

# ---------------------------------------------------------------------------
# Table S3: best single and pairwise treatment under BOTH H1299 calibrations
#   (fit including the first day, H0/PD; fit excluding the first day, Hγ/Coordination)
#   Treatment p -> (1-d)p for T = 15 days from each fit's own start, 1:1 co-culture IC.
#   Exponential growth (no carrying capacity, K = inf), as in the calibration.
#   R_N = (N_u - N_t)/N_u,  S_L = 1 - x(T),  J = (R_N + S_L)/2,  R_L = 1 - L_t/L_u (L = xN)
# ---------------------------------------------------------------------------
T, K = 15.0, np.inf      # no carrying capacity, as in the fit
TG = ["alpha", "beta", "gamma", "delta"]
SYM = {"alpha": "Rxα", "beta": "Rxβ", "gamma": "Rxγ", "delta": "Rxδ"}

def simulate(a, b, g, d, x0, n0, k=K, T=T, dt=0.01):
    """RK4 on (x, ln N) for dx/dt = x(1-x)(fL-fF)(1-N/K), dN/dt = N(x fL + (1-x) fF)(1-N/K).
    With K = inf the factor (1-N/K) is 1 (exponential growth).
    All parameters broadcast, so many doses are solved at once. Returns x(T), N(T)."""
    a, b, g, d, x = [np.array(v, float) for v in np.broadcast_arrays(a, b, g, d, x0)]
    u, lnk = np.full_like(x, np.log(n0)), np.log(k)
    def rhs(x, u):
        h = 1.0 - np.exp(u - lnk)
        fl, ff = a * x + b * (1 - x), g * x + d * (1 - x)
        return x * (1 - x) * (fl - ff) * h, (x * fl + (1 - x) * ff) * h
    for _ in range(int(round(T / dt))):
        k1 = rhs(x, u); k2 = rhs(x + dt/2*k1[0], u + dt/2*k1[1])
        k3 = rhs(x + dt/2*k2[0], u + dt/2*k2[1]); k4 = rhs(x + dt*k3[0], u + dt*k3[1])
        x = np.clip(x + dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]), 0, 1)
        u = np.minimum(u + dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]), lnk)
    return x, np.exp(u)

# ---- fitted parameters (calibrated_params_both.mat) and 1:1 initial conditions ----
FITS = {  # (alpha, beta, gamma, delta), x0 = L0/N0, N0
    ("Including first day", "untreated"): ((0.3796,  0.1232,  0.7641, 0.4475), 4.5e5/8.7e5, 8.7e5),
    ("Including first day", "+TGF-β1"):   ((0.1212,  0.1232,  0.7641, 0.4475), 4.5e5/8.7e5, 8.7e5),
    ("Excluding first day", "untreated"): ((0.3656, -0.1457,  0.0082, 0.5727), 7.6e5/1.7e6, 1.7e6),
    ("Excluding first day", "+TGF-β1"):   ((0.0811, -0.1457, -0.2403, 0.5727), 8.1e5/1.9e6, 1.9e6),
}
params = pd.DataFrame([dict(fit=f, condition=c, alpha=p[0], beta=p[1], gamma=p[2], delta=p[3],
                            x0=x0, N0=n0, x_star=(p[1]-p[3]) / ((p[1]-p[3]) + (p[2]-p[0])))
                       for (f, c), (p, x0, n0) in FITS.items()])
print("Fitted parameters and initial conditions\n", params.round(4).to_string(index=False), "\n")

# ---- dose search: singles on 101 doses, pairs on a 51 x 51 grid ----
def run(p, x0, n0, doses):
    q = {k: (1 - np.asarray(doses.get(k, 0.0))) * v for k, v in zip(TG, p)}
    return simulate(q["alpha"], q["beta"], q["gamma"], q["delta"], x0, n0)

def metrics(xt, nt, nu, Lu):
    RN = (nu - nt) / nu
    return RN, xt, 1 - xt * nt / Lu, 0.5 * (RN + 1 - xt)

d1, d2 = np.linspace(0, 1, 101), np.linspace(0, 1, 51)
results = []
for (fit, cond), (p, x0, n0) in FITS.items():
    xu, nu = run(p, x0, n0, {}); Lu = xu * nu
    print(f"{fit:20s} {cond:10s} untreated x(T) = {xu:.4f}")
    for t in TG:
        RN, xT, RL, J = metrics(*run(p, x0, n0, {t: d1}), nu, Lu); i = J.argmax()
        results.append(dict(fit=fit, condition=cond, kind="single", treatment=SYM[t], doses=f"{d1[i]:.2f}",
                            R_N=100*RN[i], x_T=xT[i], R_L=100*RL[i], J=J[i]))
    for u, v in combinations(TG, 2):
        RN, xT, RL, J = metrics(*run(p, x0, n0, {u: d2[:, None], v: d2[None, :]}), nu, Lu)
        i, j = np.unravel_index(J.argmax(), J.shape)
        results.append(dict(fit=fit, condition=cond, kind="pair", treatment=f"{SYM[u]}+{SYM[v]}",
                            doses=f"{d2[i]:.2f}, {d2[j]:.2f}", R_N=100*RN[i, j], x_T=xT[i, j],
                            R_L=100*RL[i, j], J=J[i, j]))
res = pd.DataFrame(results)
res.to_csv("TableS3_all_treatments.csv", index=False)

# ---- all treatments ranked by J ----
rnd = {"R_N": 1, "x_T": 4, "R_L": 1, "J": 3}
print("\nAll treatments, ranked by J")
print(res.sort_values(["fit", "condition", "kind", "J"], ascending=[True, False, False, False])
         .round(rnd).to_string(index=False))

# ---- Table S3 (best single, pair 1, pair 2, and Rxβ+Rxδ only if not already listed) ----
rows = []
for (fit, cond), g in res.groupby(["fit", "condition"], sort=False):
    s = g[g.kind == "single"].sort_values("J", ascending=False)
    pr = g[g.kind == "pair"].sort_values("J", ascending=False)
    picks = [("Single", s.iloc[0]), ("Pair 1", pr.iloc[0]), ("Pair 2", pr.iloc[1])]
    if "Rxβ+Rxδ" not in (pr.iloc[0].treatment, pr.iloc[1].treatment):
        picks.append(("", g[g.treatment == "Rxβ+Rxδ"].iloc[0]))
    for lab, r in picks:
        rows.append(dict(fit=fit, condition=cond, rank=lab, **r.drop(["fit", "condition", "kind"])))
S3 = pd.DataFrame(rows).round({"R_N": 1, "x_T": 3, "R_L": 1, "J": 3})
S3.to_csv("TableS3_treatment_both_fits.csv", index=False)
print(S3.to_string(index=False))

# ---- LaTeX body of Table S3 (paste between \midrule and \bottomrule) ----
TEX = {"Rxα": r"R_x\alpha", "Rxβ": r"R_x\beta", "Rxγ": r"R_x\gamma", "Rxδ": r"R_x\delta"}
tex = lambda name: "$" + "+".join(TEX[t] for t in name.split("+")) + "$"
num = lambda v: f"${v:.1f}$" if v < 0 else f"{v:.1f}"
FITLAB = {"Including first day": ("Including", "first day"), "Excluding first day": ("Excluding", "first day")}
lines, prev_fit = [], None
for (fit, cond), g in S3.groupby(["fit", "condition"], sort=False):
    if prev_fit is not None:
        lines.append(r"\midrule" if fit != prev_fit else r"\cmidrule(l){2-8}")
    for k, r in enumerate(g.itertuples()):
        first = fit != prev_fit
        fcol = FITLAB[fit][k] if first and k < 2 else ""
        ccol = ("untreated" if cond == "untreated" else r"+TGF-$\beta$1") if k == 0 else ""
        trt = f"{r.rank}: {tex(r.treatment)}" if r.rank else tex(r.treatment)
        lines.append(f"{fcol} & {ccol} & {trt} & {r.doses} & {r.R_N:.1f} & {r.x_T:.3f} & {num(r.R_L)} & {r.J:.3f} \\\\")
    prev_fit = fit
print("\n".join(lines))
