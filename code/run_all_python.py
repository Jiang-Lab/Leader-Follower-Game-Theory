"""Run every Python figure notebook in order (from the code/ folder):
    python run_all_python.py
Figures are written to ../figures. Fig. 5 and Figs. S5/S6 take several minutes each.
FigS11 reads calibrated_params.mat if the MATLAB step (run_all.m) has been run first."""
import subprocess, glob, sys
notebooks = sorted(glob.glob("Fig[1-5]_*.ipynb")) + sorted(glob.glob("FigS[1-9]_*.ipynb"),
            key=lambda s: int(s[4:].split("_")[0])) + ["FigS11_profile_likelihood_identifiability.ipynb"]
for nb in notebooks:
    print("running", nb, flush=True)
    r = subprocess.run([sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute",
                        "--inplace", "--ExecutePreprocessor.timeout=3600", nb])
    if r.returncode:
        sys.exit(f"failed: {nb}")
print("All Python figures saved in ../figures")
