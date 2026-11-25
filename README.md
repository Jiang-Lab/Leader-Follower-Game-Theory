# Growth-limited leader–follower evolutionary games reveal regime-dependent treatment trade-offs

Code and data for the manuscript by Rauf, Zhang, Robinson, Mouw, Marcus and Jiang.

A two-phenotype evolutionary game describes leader (L) and follower (F) cancer cells with payoff matrix
(α, β, γ, δ) and a shared logistic carrying capacity K:

- leader fitness f_L(x) = αx + β(1−x), follower fitness f_F(x) = γx + δ(1−x), x = L/(L+F)
- dx/dt = x(1−x)[f_L − f_F](1 − N/K), dN/dt = N[x f_L + (1−x) f_F](1 − N/K)
- interior equilibrium x* = (β−δ)/[(β−δ) + (γ−α)]
- treatment p → (1−d)p with dose d ∈ [0, 1]

The model is calibrated to H1299 NSCLC leader/follower monoculture and co-culture growth data with and
without TGF-β1, in two ways: using the full time course (days 0–4) and excluding the first day after
plating (days 1–4). The calibration uses exponential growth (no carrying capacity).

## Repository structure

```
code/      one script or notebook per figure or table, named Fig<number>_<purpose> or Table<number>_<purpose>
data/      growth data used for calibration (CSV)
figures/   output folder (filled by the scripts)
```

## Figures, tables and scripts

### Theory (Python notebooks, `code/`)

| Figure | Script | Purpose |
|---|---|---|
| Fig. 1 | `Fig1_payoff_structure_schematic.ipynb` | Payoff structure of the leader–follower game |
| Fig. 2 | `Fig2_dynamics_four_game_regimes.ipynb` | Dynamics in the four game regimes |
| Fig. 3 | `Fig3_single_treatment_tradeoffs.ipynb` | Single-treatment trade-offs between tumor reduction and leader suppression |
| Fig. 4 | `Fig4_combination_Rbeta_Rdelta_landscapes.ipynb` | R_xβ + R_xδ landscapes at T = 15 days (A) and T = 30 days (B) |
| Fig. 6 | `Fig6_optimal_treatment_game_space.py` / `.ipynb` | Best single and pairwise treatment across the (β, α) plane, for both calibrations (A–B: days 0–4; C–D: days 1–4), with the H1299 points |
| Figs. S1–S6 | `FigS1_…` to `FigS5_S6_…` notebooks | Phase lines, parameter sweeps, regime boundaries, dose responses, pairwise landscapes |

Fig. 6 reuses the saved maps in `Fig6_maps_cache.npz` (about 1 minute). Set `USE_CACHE = False`
to recompute them (about 25 minutes).

### Data calibration (`code/`)

| Figure / table | Script | Purpose |
|---|---|---|
| Fig. 5 | `Fig5_calibration_monoculture_and_coculture_fits.m` | Fits α, δ to the monocultures and β, γ to all six co-cultures; panels G–H show f_L − f_F for both calibrations. **Run first**: saves `calibrated_params.mat` |
| Table 3, Table S2 | `Table3_S2_parameters_both_calibrations.m` (Python mirror: `.ipynb` with `calibration.py`) | All four TGF-β1 cases (H₀, H_β, H_γ, H_βγ) for both calibrations: estimates, bootstrap CIs, x*, game class, ΔAICc. Saves `calibrated_params_both.mat` |
| Table S3 | `TableS3_treatment_both_fits.py` / `.ipynb` | Best single and pairwise treatments for the calibrated H1299 system under both calibrations (exponential growth, 15 days from the measured 1:1 co-culture) |
| Fig. S7 | `FigS7_four_case_fits_TGFbeta_effect.m` | Fits of the four TGF-β1 cases |
| Fig. S8 | `FigS8_model_comparison_AICc_bootstrap.m` | ΔAICc of the four cases; bootstrap distributions of Δβ, Δγ |
| Fig. S9 | `FigS9_leave_one_ratio_out_prediction.m` | Leave-one-ratio-out prediction |
| Fig. S10 | `FigS10_sensitivity_and_regime_robustness.m` | One-at-a-time sensitivity and regime robustness |
| Fig. S11 | `FigS11_profile_likelihood_identifiability.py` / `.ipynb` | Profile likelihoods (practical identifiability) |

`lf_model.py` holds the shared Python model and plotting style. Each figure is saved to `figures/`
as `.png` (300 dpi) and vector `.pdf`.

## How to run

**MATLAB** (R2020a or newer; Optimization and Statistics and Machine Learning Toolboxes), from `code/`:

```matlab
run_all          % Fig5 first, then Table3_S2 and FigS7–S10
```

The bootstrap scripts use `nBoot = 1000`; set `nBoot = 100` for a quick test. `Table3_S2` takes about 30 minutes.

**Python** (3.9 or newer), from `code/`:

```bash
pip install -r ../requirements.txt
python run_all_python.py                     # all notebooks
python TableS3_treatment_both_fits.py        # Table S3 only (about 1 minute)
```

## Data

| File | Contents |
|---|---|
| `data/monoculture_growth_normalized.csv` | Leader and follower monocultures, untreated and +TGF-β1; five replicates; days 0–4; normalized to each replicate's day-0 count (`day0_count` gives the raw value) |
| `data/coculture_counts.csv` | Leader and total cell counts in co-cultures seeded at 1:1, 1:2 and 1:5 leader:follower, untreated and +TGF-β1; one replicate per ratio; days 0–4 |

The scripts contain the same values inline, so they run without reading the CSV files.
The monoculture and 1:1 co-culture data were published in Khatib et al., *Cell Reports* 44(2), 2025;
the 1:2 and 1:5 co-culture data are new to this study.

## Methods in brief

- **Monocultures:** y(t) = y₀e^{rt} fitted to the replicate mean (r ∈ [0, 1]); replicate bootstrap, B = 1000.
- **Co-cultures:** all six series fitted jointly to x(t) and N/N₀, each weighted by its pooled SD; α, δ fixed; residual bootstrap, B = 1000. Two calibrations: days 0–4 and days 1–4.
- **Model selection:** four TGF-β1 cases compared by AICc.
- **Prediction:** leave-one-ratio-out cross-validation.
- **Identifiability:** profile likelihood.
- **Treatment of the calibrated system:** p → (1−d)p for 15 days with exponential growth; single doses on 101 values, pairs on a 51 × 51 grid.

## Citation

[Add the citation once the paper is published.]

## License

Released under the [MIT License](LICENSE).
