%% Supplementary Fig. S8: No detectable effect of TGF-beta1 on the interactions (label fig:comparison)
%
%  A    dAICc of the four co-culture fitting cases relative to the best (H0)
%  B-C  residual-bootstrap distributions of the TGF-beta1-induced changes
%       dbeta = beta_TGF - beta_U and dgamma = gamma_TGF - gamma_U under H_beta_gamma
%       (both change), with 95% CI (dashed) and zero (solid).
%
%  Cases:  H0 shared beta, gamma | Hbeta only beta differs | Hgamma only gamma
%          differs | Hbg both differ. alpha and delta fixed from the monocultures.
%  Reads alpha, delta from calibrated_params.mat (Fig6 script) if present.
%  Output: ../figures/FigS8_model_comparison.png / .pdf

clear; close all; clc; warning('off','all'); rng(1);

nBoot = 1000;                         % bootstrap samples (100 for a quick test)
tData = (0:4)';
outDir = fullfile('..', 'figures');  if ~isfolder(outDir), mkdir(outDir); end

%% ---------------- Style ----------------
baseFont = 'Helvetica';
LABEL_FS = 18;  TITLE_FS = 16;  TICK_FS = 11;  LEGEND_FS = 12;  PANEL_FS = 20;  LW = 2;
histCol = [0.45 0.45 0.45];  ciCol = [0.85 0.33 0.10];
hypNames = {'H_0', 'H_\beta', 'H_\gamma', 'H_{\beta\gamma}'};
hypK     = [2 3 3 4];
set(groot, 'DefaultFigureColor','w', 'DefaultAxesFontName',baseFont, ...
    'DefaultTextFontName',baseFont, 'DefaultAxesFontSize',TICK_FS, ...
    'DefaultAxesLabelFontSizeMultiplier',1, 'DefaultAxesTitleFontSizeMultiplier',1, ...
    'DefaultAxesLineWidth',1, 'DefaultAxesTickDir','out', 'DefaultAxesBox','off');

%% ---------------- Fixed growth rates ----------------
if isfile('calibrated_params.mat')
    cal = load('calibrated_params.mat');  alphaFix = cal.alphaFix;  deltaFix = cal.deltaFix;
else
    alphaFix = [0.3792 0.1212];  deltaFix = [0.4470 0.4470];
end

%% ---------------- Co-culture data ----------------
raw = { ...
 1, [4.5e5 7.6e5 9.0e5 1.0e6 9.4e5], [8.7e5 1.7e6 2.7e6 3.4e6 4.3e6];
 1, [3.6e5 6.2e5 5.8e5 6.5e5 6.4e5], [1.0e6 2.0e6 2.8e6 3.8e6 5.3e6];
 1, [1.8e5 2.7e5 2.6e5 3.3e5 4.0e5], [9.9e5 1.9e6 3.0e6 4.4e6 7.3e6];
 2, [4.5e5 8.1e5 7.8e5 6.6e5 5.3e5], [8.7e5 1.9e6 2.5e6 2.8e6 3.5e6];
 2, [3.6e5 5.4e5 5.5e5 5.5e5 4.3e5], [1.0e6 2.0e6 3.0e6 4.0e6 5.7e6];
 2, [1.8e5 2.6e5 2.4e5 2.7e5 2.3e5], [9.9e5 2.0e6 3.0e6 4.4e6 6.7e6]};
nSer = size(raw, 1);
for i = 1:nSer
    L = raw{i,2}(:);  N = raw{i,3}(:);
    S(i).cond = raw{i,1};  S(i).x = L ./ N;  S(i).n = N ./ N(1);   %#ok<SAGROW>
    S(i).alpha = alphaFix(S(i).cond);  S(i).delta = deltaFix(S(i).cond);
end
allX = cell2mat(arrayfun(@(s) s.x(2:end), S(:), 'UniformOutput', false));
allN = cell2mat(arrayfun(@(s) s.n(2:end), S(:), 'UniformOutput', false));
W = [std(allX, 1), std(allN, 1)];
nRes = 2 * (numel(tData) - 1) * nSer;
opts = optimoptions('lsqnonlin', 'Display','off', ...
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10, 'OptimalityTolerance',1e-10);

%% ---------------- Fit the four cases, AICc ----------------
q = zeros(4,4);  aic = zeros(1,4);
for h = 1:4
    [q(h,:), sse] = fit_case(h, S, tData, W, opts);
    kk = hypK(h) + 1;
    aic(h) = nRes*log(sse/nRes) + 2*kk + 2*kk*(kk+1)/(nRes - kk - 1);
end
dAIC = aic - min(aic);

%% ---------------- Residual bootstrap under H_beta_gamma ----------------
c = sqrt(nRes / (nRes - hypK(4)));
D = nan(nBoot, 2);
fprintf('Bootstrapping H_beta_gamma (%d samples)...\n', nBoot);
for b = 1:nBoot
    try
        p = fit_case(4, boot_resample(S, q(4,:), tData, c), tData, W, opts);
        D(b,:) = [p(3) - p(1), p(4) - p(2)];
    catch
    end
    if mod(b, 100) == 0, fprintf('  %d / %d\n', b, nBoot); end
end
D = D(~any(isnan(D), 2), :);
ciD = prctile(D, [2.5 97.5], 1);
fprintf('dAICc: %s\n', mat2str(round(dAIC, 2)));
fprintf('dbeta  95%% CI [%.3f, %.3f];  dgamma 95%% CI [%.3f, %.3f]\n', ciD(:,1), ciD(:,2));

%% ---------------- Figure ----------------
fig = figure('Units','inches', 'Position',[1 1 14 4.8]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','loose');

ax = nexttile(tl);  hold(ax,'on');  style_axes(ax, TICK_FS);
bar(ax, 1:4, dAIC, 0.6, 'FaceColor',histCol, 'EdgeColor','none');
set(ax, 'XTick',1:4, 'XTickLabel',hypNames, 'TickLabelInterpreter','tex');
xlabel(ax, 'Case', 'FontWeight','bold', 'FontSize',LABEL_FS);
ylabel(ax, '\DeltaAICc', 'FontWeight','bold', 'FontSize',LABEL_FS);
panel_label(ax, 'A', PANEL_FS);

labs = {'\Delta\beta (day^{-1})', '\Delta\gamma (day^{-1})'};
for j = 1:2
    ax = nexttile(tl);  hold(ax,'on');  style_axes(ax, TICK_FS);
    histogram(ax, D(:,j), 30, 'FaceColor',histCol, 'EdgeColor','none', ...
              'Normalization','probability', 'DisplayName','Bootstrap');
    xline(ax, 0, 'k-', 'LineWidth',LW, 'DisplayName','Zero');
    xline(ax, ciD(1,j), '--', 'Color',ciCol, 'LineWidth',LW, 'DisplayName','95% CI');
    xline(ax, ciD(2,j), '--', 'Color',ciCol, 'LineWidth',LW, 'HandleVisibility','off');
    xlabel(ax, labs{j}, 'FontWeight','bold', 'FontSize',LABEL_FS);
    if j == 1
        ylabel(ax, 'Probability', 'FontWeight','bold', 'FontSize',LABEL_FS);
        lg = legend(ax, 'Location','northwest', 'Box','off');  lg.FontSize = LEGEND_FS;
    end
    panel_label(ax, char('A' + j), PANEL_FS);
end
drawnow;
exportgraphics(fig, fullfile(outDir, 'FigS8_model_comparison.png'), 'Resolution',300);
exportgraphics(fig, fullfile(outDir, 'FigS8_model_comparison.pdf'), 'ContentType','vector');
fprintf('Saved FigS8_model_comparison.png / .pdf\n');

%% ================= Local functions =================
function q = expand_params(p, h)
    switch h
        case 1, q = [p(1) p(2) p(1) p(2)];
        case 2, q = [p(1) p(2) p(3) p(2)];
        case 3, q = [p(1) p(2) p(1) p(3)];
        otherwise, q = p(1:4);
    end
end

function [q, sse] = fit_case(h, S, t, W, opts)
    k = [2 3 3 4];  k = k(h);
    fun = @(p) case_resid(expand_params(p, h), S, t, W);
    p   = lsqnonlin(fun, 0.3*ones(1,k), -2*ones(1,k), 3*ones(1,k), opts);
    q   = expand_params(p, h);  sse = sum(fun(p).^2);
end

function r = case_resid(q, S, t, W)
    r = [];
    for i = 1:numel(S)
        [x, n] = lf_sim(q(2*S(i).cond-1 : 2*S(i).cond), S(i), t);
        r = [r; (x(2:end)-S(i).x(2:end))/W(1); (n(2:end)-S(i).n(2:end))/W(2)]; %#ok<AGROW>
    end
end

function Sb = boot_resample(S, q, t, c)
    Sb = S;
    for i = 1:numel(S)
        [x, n] = lf_sim(q(2*S(i).cond-1 : 2*S(i).cond), S(i), t);
        rx = S(i).x(2:end) - x(2:end);  rx = (rx - mean(rx)) * c;
        rn = S(i).n(2:end) - n(2:end);  rn = (rn - mean(rn)) * c;
        m = numel(rx);
        Sb(i).x = [S(i).x(1); x(2:end) + rx(randi(m, m, 1))];
        Sb(i).n = [1;         n(2:end) + rn(randi(m, m, 1))];
    end
end

function dz = rhs_2D(~, z, p)
    A = p(3)*z(1) + p(1)*(1-z(1));  B = p(2)*z(1) + p(4)*(1-z(1));
    dz = [(A - B)*z(1)*(1-z(1)); (A*z(1) + B*(1-z(1)))*z(2)];
end

function [x, n] = lf_sim(bg, Si, tout)
    p = [bg(1), bg(2), Si.alpha, Si.delta];
    [~, Z] = ode45(@(tt,z) rhs_2D(tt, z, p), tout, [Si.x(1); 1], odeset('RelTol',1e-8,'AbsTol',1e-10));
    x = Z(:,1);  n = Z(:,2);
end

function style_axes(ax, fs)
    set(ax, 'FontSize',fs, 'LineWidth',1, 'TickDir','out', 'Box','off', 'Layer','top');
    grid(ax, 'on');  ax.GridAlpha = 0.18;
end

function panel_label(ax, letter, fs)
    text(ax, -0.12, 1.04, letter, 'Units','normalized', 'FontSize',fs, 'FontWeight','bold', ...
         'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Clipping','off');
end
