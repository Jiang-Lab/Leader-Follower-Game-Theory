%% Run all MATLAB figure scripts in order (run from the code/ folder).
%  Fig6 must run first: it fits the monocultures and saves calibrated_params.mat,
%  which the supplementary scripts read so that every figure uses the same alpha and delta.
%  FigS11 (identifiability) is Python: run FigS11_profile_likelihood_identifiability.py
%  or the notebook afterwards, from this folder.
scripts = {'Fig6_calibration_monoculture_and_coculture_fits', ...
           'FigS7_four_case_fits_TGFbeta_effect', ...
           'FigS8_model_comparison_AICc_bootstrap', ...
           'FigS9_leave_one_ratio_out_prediction', ...
           'FigS10_sensitivity_and_regime_robustness'};
for k = 1:numel(scripts)
    fprintf('\n===== %s =====\n', scripts{k});
    run(scripts{k});
    close all;
end
fprintf('\nAll MATLAB figures saved in ../figures\n');
