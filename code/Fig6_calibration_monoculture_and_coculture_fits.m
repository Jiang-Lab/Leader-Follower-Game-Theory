%% Fig. 6 (label fig:fittingtgf): Calibration of the leader-follower game model
%
%  A-B  Monoculture growth of leaders (A) and followers (B), untreated (blue)
%       and +TGF-beta1 (orange), normalized to day 0. Mean +/- SD of five
%       replicates, exponential fits, 95% replicate-bootstrap bands.
%       -> alpha (leaders) and delta (followers)
%  C-D  Leader fraction x(t) in co-cultures seeded at 1:1 (magenta),
%       1:2 (green) and 1:5 (blue) leader:follower, untreated (C), +TGF (D)
%  E-F  Total cells N/N0 (log y-axis) for the same co-cultures
%  C-F  One joint fit of the game model to all six co-cultures, alpha and
%       delta fixed from A-B, beta and gamma shared by both conditions (H0),
%       95% residual-bootstrap bands.
%
%  Model:  A = alpha*x + beta*(1-x),  B = gamma*x + delta*(1-x)
%          dx/dt = (A - B) x (1-x),   dN/dt = (A x + B (1-x)) N
%
%  Output: ../figures/Fig6_calibration.png / .pdf, and calibrated_params.mat
%  (read by the FigS7-S11 scripts). RUN THIS SCRIPT FIRST.
%  Needs:  Optimization Toolbox (lsqcurvefit, lsqnonlin), Statistics Toolbox
%          (prctile, randsample), MATLAB R2020a or newer.

clear; close all; clc; warning('off','all'); rng(1);   % reproducible

%% ---------------- Settings ----------------
nBoot = 1000;                        % bootstrap samples (lower for a quick test)
tData = (0:4)';                      % days
tGrid = linspace(0, 4, 101)';

%% ---------------- Style (PNAS-like, as in the monoculture figures) ----------------
baseFont       = 'Helvetica';
LABEL_FS       = 18;                 % axis labels (bold)
TITLE_FS       = 16;                 % panel titles (bold)
TICK_FS        = 11;                 % tick labels
LEGEND_FS      = 12;                 % legends
PANEL_FS       = 20;                 % panel letters (bold)
FIT_LW         = 2;                  % fitted curves
DATA_LW        = 1.5;                % markers / error bars
MARKER_SIZE    = 6;
BAND_ALPHA     = 0.12;

colorOrder = [0.95 0.30 0.95;        % magenta
              0.18 0.60 0.25;        % green
              0.00 0.45 0.74;        % blue
              0.85 0.33 0.10;        % orange (red-orange)
              0.95 0.69 0.15];       % yellow-orange
condCol  = colorOrder([3 4], :);     % untreated blue, +TGF-beta1 orange
ratioCol = colorOrder(1:3, :);       % 1:1 magenta, 1:2 green, 1:5 blue
ratioMk  = {'o', 's', '^'};
ratioLabels = {'1:1 L:F', '1:2 L:F', '1:5 L:F'};
condLabels  = {'Untreated', '+TGF-\beta1'};

set(groot, 'DefaultFigureColor','w', 'DefaultAxesFontName',baseFont, ...
    'DefaultTextFontName',baseFont, 'DefaultAxesFontSize',TICK_FS, ...
    'DefaultAxesLabelFontSizeMultiplier',1, 'DefaultAxesTitleFontSizeMultiplier',1, ...
    'DefaultAxesLineWidth',1, 'DefaultLineLineWidth',FIT_LW, ...
    'DefaultAxesTickDir','out', 'DefaultAxesBox','off', ...
    'DefaultAxesGridLineStyle','-', 'DefaultAxesGridAlpha',0.18);

%% ================= A-B: monoculture data and fits =================
% Normalized N/N0; rows = days 0-4, columns = 5 replicates
% Leaders   raw day-0: 49000, 43800, 46400, 46400, 8.9E+05
% Followers raw day-0: 52000, 68500, 54200, 54200, 1.2E+06
mono{1,1} = [1.0000 1.0000 1.0000 1.0000 1.0000;     % leaders, untreated
             2.3510 2.4155 1.1853 1.2026 2.1348;
             3.0204 3.2785 2.8707 3.0000 2.8090;
             4.5306 5.9361 3.6293 4.2716 3.5955;
             6.8163 6.4840 5.5603 5.9052 4.4944];
mono{1,2} = [1.0000 1.0000 1.0000 1.0000 1.0000;     % leaders, +TGF-beta1
             0.8143 0.7173 1.7026 1.7974 1.4607;
             1.4051 0.8342 2.2586 2.3147 1.6854;
             0.6257 0.9541 2.0560 2.3534 1.7978;
             1.0212 1.0425 2.7026 2.4828 1.9101];
mono{2,1} = [1.0000 1.0000 1.0000 1.0000 1.0000;     % followers, untreated
             1.2577 1.5942 1.5055 1.5240 1.8333;
             1.7298 1.0526 1.2362 1.1107 2.7500;
             3.5769 2.8321 3.8376 3.8376 4.1667;
             6.6000 4.9861 6.7159 6.3100 6.6667];
mono{2,2} = [1.0000 1.0000 1.0000 1.0000 1.0000;     % followers, +TGF-beta1
             1.7077 0.9530 0.9041 0.9041 1.7500;
             1.3154 0.9372 1.2362 1.2694 2.6667;
             4.0192 3.0511 3.2657 3.2030 3.8333;
             4.7731 3.5365 6.2362 6.4945 6.5000];
monoNames = {'Leaders', 'Followers'};

optsMono = optimoptions('lsqcurvefit', 'Display','off', ...
    'FunctionTolerance',1e-12, 'StepTolerance',1e-12, 'OptimalityTolerance',1e-12);

tic;
fprintf('Fitting monocultures (%d bootstrap samples)...\n', nBoot);
for ct = 1:2
    for cnd = 1:2
        M(ct,cnd) = mono_fit(mono{ct,cnd}, tData, tGrid, nBoot, optsMono); %#ok<SAGROW>
    end
end
alphaFix = [M(1,1).rMean, M(1,2).rMean];   % leaders:   untreated, +TGF
deltaFix = [M(2,1).rMean, M(2,1).rMean];   % followers: one shared value (untreated)

%% ================= C-F: co-culture data and joint fit (H0) =================
% condition (1 untreated, 2 +TGF), ratio index (1 = 1:1, 2 = 1:2, 3 = 1:5 L:F),
% leader counts, total counts
raw = { ...
 1, 1, [4.5e5 7.6e5 9.0e5 1.0e6 9.4e5], [8.7e5 1.7e6 2.7e6 3.4e6 4.3e6];
 1, 2, [3.6e5 6.2e5 5.8e5 6.5e5 6.4e5], [1.0e6 2.0e6 2.8e6 3.8e6 5.3e6];
 1, 3, [1.8e5 2.7e5 2.6e5 3.3e5 4.0e5], [9.9e5 1.9e6 3.0e6 4.4e6 7.3e6];
 2, 1, [4.5e5 8.1e5 7.8e5 6.6e5 5.3e5], [8.7e5 1.9e6 2.5e6 2.8e6 3.5e6];
 2, 2, [3.6e5 5.4e5 5.5e5 5.5e5 4.3e5], [1.0e6 2.0e6 3.0e6 4.0e6 5.7e6];
 2, 3, [1.8e5 2.6e5 2.4e5 2.7e5 2.3e5], [9.9e5 2.0e6 3.0e6 4.4e6 6.7e6]};
nSer = size(raw, 1);
for i = 1:nSer
    L = raw{i,3}(:);  N = raw{i,4}(:);
    S(i).cond  = raw{i,1};          %#ok<SAGROW>
    S(i).ratio = raw{i,2};
    S(i).x     = L ./ N;             % leader fraction
    S(i).n     = N ./ N(1);          % N/N0
    S(i).alpha = alphaFix(S(i).cond);
    S(i).delta = deltaFix(S(i).cond);
end

% Residual scaling: pooled SD of x and of N/N0 over days 1-4, all series
allX = cell2mat(arrayfun(@(s) s.x(2:end), S(:), 'UniformOutput', false));
allN = cell2mat(arrayfun(@(s) s.n(2:end), S(:), 'UniformOutput', false));
W = [std(allX, 1), std(allN, 1)];
nRes = 2 * (numel(tData) - 1) * nSer;              % 48 fitted points

opts = optimoptions('lsqnonlin', 'Display','off', ...
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10, 'OptimalityTolerance',1e-10);

bgHat = fit_shared(S, tData, W, opts);             % [beta, gamma]

fprintf('Bootstrapping the co-culture fit (%d samples)...\n', nBoot);
c  = sqrt(nRes / (nRes - 2));                     % small-sample inflation
PB = nan(nBoot, 2);
for b = 1:nBoot
    try
        PB(b,:) = fit_shared(boot_resample(S, bgHat, tData, c), tData, W, opts);
    catch
    end
    if mod(b, 100) == 0, fprintf('  %d / %d\n', b, nBoot); end
end
PB  = PB(~any(isnan(PB), 2), :);
ciB = prctile(PB, [2.5 97.5], 1);

% Fitted curves and bands for every co-culture
for i = 1:nSer
    [xx, nn] = lf_sim(bgHat, S(i), tGrid);
    S(i).xFit = xx;  S(i).nFit = nn;
    Xb = zeros(size(PB,1), numel(tGrid));  Nb = Xb;
    for b = 1:size(PB,1)
        [xx, nn] = lf_sim(PB(b,:), S(i), tGrid);
        Xb(b,:) = xx';  Nb(b,:) = nn';
    end
    S(i).xCI = prctile(Xb, [2.5 97.5], 1);
    S(i).nCI = prctile(Nb, [2.5 97.5], 1);
end

%% ---------------- Numbers for the caption ----------------
fprintf('\nMonoculture growth rates (bootstrap mean, 95%% CI, 1/day):\n');
for ct = 1:2
    for cnd = 1:2
        fprintf('  %-9s %-12s %.4f  [%.4f, %.4f]\n', monoNames{ct}, ...
            strrep(condLabels{cnd},'\',''), M(ct,cnd).rMean, M(ct,cnd).rCI(1), M(ct,cnd).rCI(2));
    end
end
fprintf('Used in co-culture fit: alpha = %.4f / %.4f, delta = %.4f (shared)\n', ...
        alphaFix(1), alphaFix(2), deltaFix(1));
fprintf('Shared interactions (H0): beta  = %.4f [%.4f, %.4f]\n', bgHat(1), ciB(1,1), ciB(2,1));
fprintf('                          gamma = %.4f [%.4f, %.4f]\n', bgHat(2), ciB(1,2), ciB(2,2));

% Save the calibrated parameters for the supplementary-figure scripts (FigS7-S11)
rCI_alpha = [M(1,1).rCI; M(1,2).rCI];      % rows: untreated, +TGF
rCI_delta = M(2,1).rCI;
save('calibrated_params.mat', 'alphaFix', 'deltaFix', 'bgHat', 'ciB', 'rCI_alpha', 'rCI_delta');
fprintf('Saved calibrated_params.mat\n');

%% ================= Figure =================
fig = figure('Units','inches', 'Position',[1 1 10 11.5]);
tl  = tiledlayout(fig, 3, 2, 'TileSpacing','compact', 'Padding','compact');
letters = 'ABCDEF';

% --- A-B: monocultures ---
for ct = 1:2
    ax = nexttile(tl, ct);  hold(ax, 'on');  style_axes(ax, TICK_FS);
    h = gobjects(1, 2);
    for cnd = 1:2
        col = condCol(cnd,:);  R = M(ct,cnd);
        add_band(ax, tGrid, R.band(1,:), R.band(2,:), col, BAND_ALPHA);
        h(cnd) = plot(ax, tGrid, R.fit, '-', 'Color',col, 'LineWidth',FIT_LW);
        errorbar(ax, tData, mean(mono{ct,cnd},2), std(mono{ct,cnd},0,2), 'o', ...
            'Color',col, 'MarkerFaceColor',col, 'MarkerSize',MARKER_SIZE, ...
            'LineWidth',DATA_LW, 'CapSize',6, 'HandleVisibility','off');
    end
    title(ax, [monoNames{ct} ' (Monoculture)'], 'FontWeight','bold', 'FontSize',TITLE_FS);
    ylabel(ax, [monoNames{ct}(1:end-1) ' Cells (N/N_0)'], 'FontWeight','bold', 'FontSize',LABEL_FS);
    xlabel(ax, 'Time (days)', 'FontWeight','bold', 'FontSize',LABEL_FS);
    if ct == 1
        lg = legend(ax, h, condLabels, 'Location','northwest', 'Box','off');
        lg.FontSize = LEGEND_FS;  lg.ItemTokenSize = [16 10];
    end
    panel_label(ax, letters(ct), PANEL_FS);
end

% --- C-F: co-cultures (left column untreated, right column +TGF-beta1) ---
for row = 1:2                                  % 1 = leader fraction, 2 = N/N0
    for cnd = 1:2
        tileNo = 2 + (row-1)*2 + cnd;
        ax = nexttile(tl, tileNo);  hold(ax, 'on');  style_axes(ax, TICK_FS);
        h = gobjects(1, 3);
        for i = find([S.cond] == cnd)
            rr = S(i).ratio;  col = ratioCol(rr,:);
            if row == 1
                yF = S(i).xFit;  yC = S(i).xCI;  yD = S(i).x;
            else
                yF = S(i).nFit;  yC = S(i).nCI;  yD = S(i).n;
            end
            add_band(ax, tGrid, yC(1,:), yC(2,:), col, BAND_ALPHA);
            h(rr) = plot(ax, tGrid, yF, '-', 'Color',col, 'LineWidth',FIT_LW);
            plot(ax, tData, yD, ratioMk{rr}, 'Color',col, 'MarkerFaceColor',col, ...
                 'MarkerSize',MARKER_SIZE, 'LineWidth',DATA_LW, 'HandleVisibility','off');
        end
        if row == 1
            title(ax, ['Co-culture, ' condLabels{cnd}], 'FontWeight','bold', 'FontSize',TITLE_FS);
            ylabel(ax, 'Leader Fraction x(t)', 'FontWeight','bold', 'FontSize',LABEL_FS);
            ylim(ax, [0 0.6]);
        else
            set(ax, 'YScale','log');
            ylabel(ax, 'Total Cells (N/N_0)', 'FontWeight','bold', 'FontSize',LABEL_FS);
            xlabel(ax, 'Time (days)', 'FontWeight','bold', 'FontSize',LABEL_FS);
        end
        if tileNo == 3
            lg = legend(ax, h, ratioLabels, 'Location','northeast', 'Box','off');
            lg.FontSize = LEGEND_FS;  lg.ItemTokenSize = [16 10];
        end
        panel_label(ax, letters(tileNo), PANEL_FS);
    end
end
drawnow;                                    % show the finished figure before saving
fprintf('Saving figure...\n');
outDir = fullfile('..', 'figures');  if ~isfolder(outDir), mkdir(outDir); end
exportgraphics(fig, fullfile(outDir, 'Fig6_calibration.png'), 'Resolution',300);
exportgraphics(fig, fullfile(outDir, 'Fig6_calibration.pdf'), 'ContentType','vector');
fprintf('\nSaved Fig6_calibration.png and .pdf (total time %.1f min)\n', toc/60);

%% ================= Local functions =================
function R = mono_fit(T, t, tGrid, nBoot, opts)
% Exponential fit y = y0*exp(r t) to the replicate mean (exact solution of
% y' = r y), bounds lb = [0, min(T)], ub = [1, max(T)], replicate bootstrap.
    model = @(p,tt) p(2) .* exp(p(1) .* tt);
    lb = [0.0, min(T(:))];  ub = [1.0, max(T(:))];  p0 = [0.5, T(1,1)];
    nRep = size(T, 2);  P = nan(nBoot, 2);
    for b = 1:nBoot
        idx = randsample(nRep, nRep, true);
        try
            P(b,:) = lsqcurvefit(model, p0, t, mean(T(:,idx), 2), lb, ub, opts);
        catch
        end
    end
    P = P(~any(isnan(P), 2), :);
    Y = P(:,2) .* exp(P(:,1) .* tGrid');
    R.rMean = mean(P(:,1));
    R.rCI   = prctile(P(:,1), [2.5 97.5]);
    R.fit   = mean(Y, 1);
    R.band  = prctile(Y, [2.5 97.5], 1);
end

function dz = rhs_2D(~, z, p)
    beta = p(1);  gamma = p(2);  alpha = p(3);  delta = p(4);
    x = z(1);  N = z(2);
    A = alpha*x + beta*(1-x);
    B = gamma*x + delta*(1-x);
    dz = [(A - B)*x*(1-x);
          (A*x + B*(1-x))*N];
end

function [x, n] = lf_sim(bg, Si, tout)
    p = [bg(1), bg(2), Si.alpha, Si.delta];
    [~, Z] = ode45(@(tt,z) rhs_2D(tt, z, p), tout, [Si.x(1); 1], ...
                   odeset('RelTol',1e-8, 'AbsTol',1e-10));
    x = Z(:,1);  n = Z(:,2);
end

function r = lf_resid(bg, S, t, W)
    r = [];
    for i = 1:numel(S)
        [x, n] = lf_sim(bg, S(i), t);
        r = [r; (x(2:end) - S(i).x(2:end)) / W(1);
                (n(2:end) - S(i).n(2:end)) / W(2)]; %#ok<AGROW>
    end
end

function bg = fit_shared(S, t, W, opts)
% Joint fit of beta and gamma shared by all six co-cultures (case H0)
    bg = lsqnonlin(@(p) lf_resid(p, S, t, W), [0.3 0.3], [-2 -2], [3 3], opts);
end

function Sb = boot_resample(S, bg, t, c)
% Residual bootstrap for single-replicate series
    Sb = S;
    for i = 1:numel(S)
        [x, n] = lf_sim(bg, S(i), t);
        rx = S(i).x(2:end) - x(2:end);   rx = (rx - mean(rx)) * c;
        rn = S(i).n(2:end) - n(2:end);   rn = (rn - mean(rn)) * c;
        m  = numel(rx);
        Sb(i).x = [S(i).x(1); x(2:end) + rx(randi(m, m, 1))];
        Sb(i).n = [1;         n(2:end) + rn(randi(m, m, 1))];
    end
end

function style_axes(ax, fs)
    set(ax, 'FontSize',fs, 'LineWidth',1, 'TickDir','out', 'Box','off', 'Layer','top');
    grid(ax, 'on');  ax.GridAlpha = 0.18;
end

function add_band(ax, t, lo, hi, col, a)
    fill(ax, [t; flipud(t)], [lo(:); flipud(hi(:))], col, ...
         'FaceAlpha',a, 'EdgeColor','none', 'HandleVisibility','off');
end

function panel_label(ax, letter, fs)
    text(ax, -0.12, 1.04, letter, 'Units','normalized', 'FontSize',fs, ...
         'FontWeight','bold', 'HorizontalAlignment','left', ...
         'VerticalAlignment','bottom', 'Clipping','off');
end
