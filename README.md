# Growth-limited leader–follower evolutionary games reveal regime-dependent treatment trade-offs

Code and data for the manuscript by Rauf, Zhang, Robinson, Mouw, Marcus and Jiang.

A two-phenotype evolutionary game describes leader (L) and follower (F) cancer cells with payoff matrix
(α, β, γ, δ) and a shared logistic carrying capacity K:

- leader fitness  f_L(x) = αx + β(1−x),  follower fitness  f_F(x) = γx + δ(1−x),  x = L/(L+F)
- dx/dt = x(1−x)[f_L − f_F](1 − N/K),  dN/dt = N[x f_L + (1−x) f_F](1 − N/K)

The model is calibrated to H1299 NSCLC leader/follower monoculture and co-culture growth data with and without TGF-β1.

## Repository structure

```
code/      figure scripts, one per figure, named Fig<number>_<purpose>
data/      growth data used for calibration (CSV)
figures/   output folder (created and filled by the scripts)
theory/    scripts for the theoretical figures (Figs. 1–5, S1–S6)
```

## Figures and scripts

### Data-calibration figures (`code/`)

| Figure | LaTeX label | Script | Purpose | Output |
|---|---|---|---|---|
| Fig. 6 | `fig:fittingtgf` | `Fig6_calibration_monoculture_and_coculture_fits.m` | Fits α, δ to monocultures (A–B) and the shared β, γ to all six co-cultures (C–F); bootstrap 95% bands. **Run first**: saves `calibrated_params.mat` | `Fig6_calibration` |
| Fig. S7 | `fig:hypfits` | `FigS7_four_case_fits_TGFbeta_effect.m` | Fits of the four cases for the TGF-β1 effect (H₀, H_β, H_γ, H_βγ) overlaid on the data | `FigS7_four_case_fits` |
| Fig. S8 | `fig:comparison` | `FigS8_model_comparison_AICc_bootstrap.m` | ΔAICc of the four cases; bootstrap distributions of Δβ, Δγ under H_βγ | `FigS8_model_comparison` |
| Fig. S9 | `fig:loo` | `FigS9_leave_one_ratio_out_prediction.m` | Leave-one-ratio-out prediction: fit two seeding ratios, predict the third | `FigS9_leave_one_ratio_out` |
| Fig. S10 | `fig:sens1` | `FigS10_sensitivity_and_regime_robustness.m` | ±10% one-at-a-time sensitivity (A–D); regime robustness in the (β−δ, γ−α) plane with ±25% Monte Carlo (E–F). Set `K_cells` at the top | `FigS10_sensitivity_K1e10` |
| Fig. S11 | `fig:ide` | `FigS11_profile_likelihood_identifiability.py` / `.ipynb` | Profile likelihoods of α, δ, β, γ, the joint (β, γ) region, and Δβ, Δγ under the alternative cases | `FigS11_identifiability` |

Each output is saved to `figures/` as both `.png` (300 dpi) and vector `.pdf`.

### Theory figures (`theory/`)

See `theory/README.md` for Figs. 1–5 and S1–S6.

## How to run

**MATLAB** (R2020a or newer; Optimization and Statistics and Machine Learning Toolboxes), from the `code/` folder:

```matlab
run_all          % runs Fig6 first, then FigS7–S10
```

or run the scripts one by one, starting with `Fig6_calibration_monoculture_and_coculture_fits.m`.
The bootstrap scripts use `nBoot = 1000` (Fig6 and FigS8: a few minutes; FigS9: about 3× longer). Set `nBoot = 100` for a quick test.

**Python** (3.9 or newer), from the `code/` folder, after the MATLAB step:

```bash
pip install -r ../requirements.txt
python FigS11_profile_likelihood_identifiability.py      # about 35 s
```

or open `FigS11_profile_likelihood_identifiability.ipynb` and run all cells.
If `calibrated_params.mat` is present it is used for α and δ; otherwise the values typed in the script are used.

## Data

| File | Contents |
|---|---|
| `data/monoculture_growth_normalized.csv` | Leader and follower monocultures, untreated and +TGF-β1; five replicates; days 0–4; cell number normalized to each replicate's day-0 count (`day0_count` gives the raw value) |
| `data/coculture_counts.csv` | Leader and total cell counts in co-cultures seeded at 1:1, 1:2 and 1:5 leader:follower, untreated and +TGF-β1; one replicate per ratio; days 0–4 |

The scripts contain the same values inline, so they run without reading the CSV files.
The monoculture and 1:1 co-culture data were published in Khatib et al., *Cell Reports* 44(2), 2025; the 1:2 and 1:5 co-culture data are new to this study.

## Methods in brief

- **Monocultures:** y(t) = y₀e^{rt} fitted to the replicate mean with `lsqcurvefit` (r ∈ [0, 1], y₀ within the data range); replicate bootstrap, B = 1000. One δ (untreated) is used for both conditions.
- **Co-cultures:** all six series fitted jointly to x(t) and N/N₀, each divided by its pooled SD, with `lsqnonlin` and `ode45`; α, δ fixed; β, γ ∈ [−2, 3]. Residual bootstrap, B = 1000.
- **Model selection:** four cases compared by AICc (n = 48).
- **Prediction:** leave-one-ratio-out cross-validation.
- **Identifiability:** profile likelihood with Gaussian errors and the variance profiled out.

## Citation

[Add the citation once the paper is published.]

## License

[Choose a license, e.g. MIT, and add it as `LICENSE`.]
