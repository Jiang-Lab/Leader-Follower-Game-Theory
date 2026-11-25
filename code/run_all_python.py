"""Run every Python notebook in order (from the code/ folder):
    python run_all_python.py
Figures are written to ../figures. Fig. 4, Fig. 6 (without the saved maps) and Figs. S5/S6
take several minutes each. FigS11 reads calibrated_params.mat if the MATLAB step
(run_all.m) has been run first; otherwise it uses the values typed in the script.
Table3_S2 (Python mirror of the MATLAB table script) and TableS3 print the table values."""
import subprocess, sys

notebooks = [
    "Fig1_payoff_structure_schematic.ipynb",
    "Fig2_dynamics_four_game_regimes.ipynb",
    "Fig3_single_treatment_tradeoffs.ipynb",
    "Fig4_combination_Rbeta_Rdelta_landscapes.ipynb",
    "Fig6_optimal_treatment_game_space.ipynb",
    "FigS1_phase_lines_growth_coupled.ipynb",
    "FigS2_parameter_sweeps_regime_map.ipynb",
    "FigS3_dynamics_on_regime_boundaries.ipynb",
    "FigS4_single_treatment_dose_response.ipynb",
    "FigS5_S6_pairwise_combination_landscapes.ipynb",
    "FigS11_profile_likelihood_identifiability.ipynb",
    "Table3_S2_parameters_both_calibrations.ipynb",
    "TableS3_treatment_both_fits.ipynb",
]
for nb in notebooks:
    print("running", nb, flush=True)
    r = subprocess.run([sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute",
                        "--inplace", "--ExecutePreprocessor.timeout=7200", nb])
    if r.returncode:
        sys.exit(f"failed: {nb}")
print("All Python figures saved in ../figures")
