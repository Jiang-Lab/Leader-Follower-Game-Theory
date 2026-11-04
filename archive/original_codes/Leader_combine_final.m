%% ==========================================================
%  ONE-PANEL overlay of TWO datasets:
%  exponential fit + bootstrap 95% band + RMSE/R2
% ==========================================================
format long; warning('off','all'); rng(1);

%% ---------- PNAS-like style ----------
baseFont = 'Helvetica';
PNAS_LABEL_FONTSIZE  = 18;
PNAS_TICK_FONTSIZE   = 11;
PNAS_LEGEND_FONTSIZE = 12;

colorOrder = [0.95 0.30 0.95;
              0.18 0.60 0.25;
              0.00 0.45 0.74;
              0.85 0.33 0.10;
              0.95 0.69 0.15];

set(groot, ...
  'DefaultFigureColor','w', ...
  'DefaultAxesFontName',baseFont, ...
  'DefaultTextFontName',baseFont, ...
  'DefaultAxesFontSize',PNAS_TICK_FONTSIZE, ...
  'DefaultAxesLabelFontSizeMultiplier',1, ...
  'DefaultAxesTitleFontSizeMultiplier',1, ...
  'DefaultAxesLineWidth',1, ...
  'DefaultLineLineWidth',2, ...
  'DefaultAxesTickDir','out', ...
  'DefaultAxesBox','off', ...
  'DefaultAxesGridLineStyle','-', ...
  'DefaultAxesGridAlpha',0.18, ...
  'DefaultAxesColorOrder',colorOrder);

%% =========================
% Shared time
% =========================
tdata = [0; 1; 2; 3; 4];
tgrid = linspace(min(tdata), max(tdata), 200);

%% =========================
% DATASET 1
% =========================
Tdata1 = [ ...
    1,          1,          1,          1,          1.0; ...
    2.351020408,2.415525114,1.185344828,1.202586207,2.1; ...
    3.020408163,3.278538813,2.870689655,3.0,        2.8; ...
    4.530612245,5.936073059,3.629310345,4.271551724,3.6; ...
    6.816326531,6.484018265,5.560344828,5.905172414,4.5 ...
];

st1 = [0; 0.595122; 0.147925135; 0.842579264; 0.488956916];

%% =========================
% DATASET 2
% =========================
Tdata2 = [ ...
    1,          1,          1,          1,          1.0; ...
    1.493877551,1.630136986,1.702586207,1.797413793,1.46; ...
    1.983673469,1.885844749,2.258620690,2.314655172,1.69; ...
    2.142857143,2.534246575,2.056034483,2.353448276,1.80; ...
    1.334693878,1.511415525,2.702586207,2.482758621,1.91 ...
];

st2 = [0; 0.126209035; 0.234470262; 0.252343355; 0.532057738];

%% =========================
% Optimizer options
% =========================
options = optimoptions('lsqcurvefit', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12);

nBoot = 1000;

%% =========================
% Fit + bootstrap each dataset separately
% =========================
res1 = boot_exp_fit_blockcols(tdata, tgrid, Tdata1, st1, nBoot, options);
res2 = boot_exp_fit_blockcols(tdata, tgrid, Tdata2, st2, nBoot, options);

%% =========================
% Print results
% =========================
disp('================ Dataset 1 ================')
fprintf('Objective-function best fit:\n');
fprintf('Best a  = %.12f\n', res1.xBest(1));
fprintf('Best y0 = %.12f\n', res1.xBest(2));
fprintf('RMSE       = %.6f\n', res1.RMSE);
fprintf('R^2        = %.6f\n', res1.R2);
fprintf('log RMSE   = %.6f\n', res1.RMSE_log);
fprintf('log R^2    = %.6f\n\n', res1.R2_log);

fprintf('Bootstrap summary:\n');
fprintf('Mean a  = %.12f | 95%% CI [%.12f, %.12f]\n', ...
    res1.paramMean(1), res1.paramCI(1,1), res1.paramCI(2,1));
fprintf('Mean y0 = %.12f | 95%% CI [%.12f, %.12f]\n\n', ...
    res1.paramMean(2), res1.paramCI(1,2), res1.paramCI(2,2));

disp('================ Dataset 2 ================')
fprintf('Objective-function best fit:\n');
fprintf('Best a  = %.12f\n', res2.xBest(1));
fprintf('Best y0 = %.12f\n', res2.xBest(2));
fprintf('RMSE       = %.6f\n', res2.RMSE);
fprintf('R^2        = %.6f\n', res2.R2);
fprintf('log RMSE   = %.6f\n', res2.RMSE_log);
fprintf('log R^2    = %.6f\n\n', res2.R2_log);

fprintf('Bootstrap summary:\n');
fprintf('Mean a  = %.12f | 95%% CI [%.12f, %.12f]\n', ...
    res2.paramMean(1), res2.paramCI(1,1), res2.paramCI(2,1));
fprintf('Mean y0 = %.12f | 95%% CI [%.12f, %.12f]\n\n', ...
    res2.paramMean(2), res2.paramCI(1,2), res2.paramCI(2,2));

%% ==========================================================
% Plot BOTH in ONE PANEL
% ==========================================================
fig = figure('Units','inches','Position',[1 1 5.6 4.2]); hold on; box on;

col1 = [0.00 0.45 0.74];
col2 = [0.85 0.33 0.10];

% ---- Dataset 1 ----
errorbar(tdata, mean(Tdata1,2), st1(:), '*', ...
    'Color', col1, ...
    'LineWidth',1.5, ...
    'MarkerSize',6, ...
    'CapSize',6, ...
    'DisplayName','Data 1 mean \pm SD');

fill([tgrid, fliplr(tgrid)], ...
     [res1.predCI(1,:), fliplr(res1.predCI(2,:))], ...
     col1, ...
     'FaceAlpha',0.10, ...
     'EdgeColor','none', ...
     'HandleVisibility','off');

plot(tgrid, res1.predBest, '-', ...
    'Color', col1, ...
    'LineWidth',2.2, ...
    'DisplayName','Fit 1');

% ---- Dataset 2 ----
errorbar(tdata, mean(Tdata2,2), st2(:), 'o', ...
    'Color', col2, ...
    'LineWidth',1.5, ...
    'MarkerSize',4.5, ...
    'CapSize',6, ...
    'DisplayName','Data 2 mean \pm SD');

fill([tgrid, fliplr(tgrid)], ...
     [res2.predCI(1,:), fliplr(res2.predCI(2,:))], ...
     col2, ...
     'FaceAlpha',0.10, ...
     'EdgeColor','none', ...
     'HandleVisibility','off');

plot(tgrid, res2.predBest, '--', ...
    'Color', col2, ...
    'LineWidth',2.2, ...
    'DisplayName','Fit 2');

xlabel('Time (days)', 'FontWeight','bold', 'FontSize',PNAS_LABEL_FONTSIZE);
ylabel('Follower Cells', 'FontWeight','bold', 'FontSize',PNAS_LABEL_FONTSIZE);

ax = gca;
grid(ax,'on');
ax.Layer = 'top';

text(ax, -0.12, 1.04, 'A', ...
    'Units','normalized', ...
    'FontSize',20, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'Clipping','off');

lg = legend('Location','eastoutside','Box','off');
lg.FontSize = PNAS_LEGEND_FONTSIZE;
lg.ItemTokenSize = [16 10];

exportgraphics(fig, 'overlay_two_datasets_onepanel.png', 'Resolution', 300);
exportgraphics(fig, 'overlay_two_datasets_onepanel.pdf', 'ContentType','vector');

%% ========================= LOCAL FUNCTIONS =========================

function dydt = g(~, y, x)
    dydt = x(1) * y;
end

function y = cost_function(x, tvec)
    [~, Y] = ode23s(@(t,y) g(t,y,x), tvec, x(2));
    y = Y(:,1);
end

function res = boot_exp_fit_blockcols(tdata, tgrid, Tdata, st, nBoot, options)
    %#ok<INUSD>

    [~, nRep] = size(Tdata);
    Ymean = mean(Tdata, 2);

    lb = [0.0, min(Tdata(:))];
    ub = [1.0, max(Tdata(:))];
    x0 = [0.5, Tdata(1,1)];

    % Objective-function best fit
    xBest = lsqcurvefit( ...
        @(x,t) cost_function(x,t), ...
        x0, ...
        tdata, ...
        Ymean, ...
        lb, ...
        ub, ...
        options);

    res.xBest = xBest;

    YfitBest = cost_function(xBest, tdata);

    [res.R2, res.RMSE] = goodness_of_fit(Ymean, YfitBest);
    [res.R2_log, res.RMSE_log] = goodness_of_fit(log(Ymean), log(YfitBest));

    % Bootstrap resampling over replicate columns
    paramEsts = nan(nBoot, numel(x0));

    for b = 1:nBoot
        idx = randsample(nRep, nRep, true);
        Tb  = Tdata(:, idx);
        syntheticData = mean(Tb, 2);

        try
            xBoot = lsqcurvefit( ...
                @(x,t) cost_function(x,t), ...
                x0, ...
                tdata, ...
                syntheticData, ...
                lb, ...
                ub, ...
                options);

            paramEsts(b,:) = xBoot(:).';

        catch
            paramEsts(b,:) = NaN;
        end
    end

    paramEsts = paramEsts(~any(isnan(paramEsts),2), :);

    res.paramMean = mean(paramEsts,1);
    res.paramCI   = squeeze(prctile(paramEsts,[2.5 97.5],1));

    % Bootstrap prediction band
    yBoot = zeros(size(paramEsts,1), numel(tgrid));

    for b = 1:size(paramEsts,1)
        yBoot(b,:) = cost_function(paramEsts(b,:), tgrid);
    end

    res.predMean = mean(yBoot,1);
    res.predCI   = squeeze(prctile(yBoot,[2.5 97.5],1));

    % Best-fit curve
    res.predBest = cost_function(xBest, tgrid);
end

function [R2, RMSE] = goodness_of_fit(yobs, yfit)

    yobs = yobs(:);
    yfit = yfit(:);

    RMSE = sqrt(mean((yobs - yfit).^2));

    SSres = sum((yobs - yfit).^2);
    SStot = sum((yobs - mean(yobs)).^2);

    R2 = 1 - SSres/SStot;
end