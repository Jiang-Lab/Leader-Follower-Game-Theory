function run_all()
%% Run all MATLAB scripts in order (run from the code/ folder).
%  Fig5 must run first: it fits the monocultures and saves calibrated_params.mat,
%  which the supplementary scripts read so that every figure uses the same alpha and delta.
%  Table3_S2 prints every value of Table 3 and Supplementary Table S2 (both calibrations,
%  all four cases; about 30 min with nBoot = 1000).
%  FigS11 (identifiability) is Python: run FigS11_profile_likelihood_identifiability.py
%  or the notebook afterwards, from this folder.
%  Each script starts with "clear", so it is run in the base workspace; the list below
%  lives in this function's workspace and is not cleared.
scripts = {'Fig5_calibration_monoculture_and_coculture_fits', ...
           'Table3_S2_parameters_both_calibrations', ...
           'FigS7_four_case_fits_TGFbeta_effect', ...
           'FigS8_model_comparison_AICc_bootstrap', ...
           'FigS9_leave_one_ratio_out_prediction', ...
           'FigS10_sensitivity_and_regime_robustness'};
for k = 1:numel(scripts)
    fprintf('\n===== %s =====\n', scripts{k});
    evalin('base', sprintf('run(''%s'')', scripts{k}));
    close all;
end
fprintf('\nAll MATLAB figures saved in ../figures\n');
end
