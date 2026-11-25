%% Supplementary Fig. S9: Leave-one-ratio-out prediction (label fig:loo)
%
%  The co-culture model (case H0: beta, gamma shared by both conditions) is
%  fitted to two seeding ratios and used, without refitting, to predict the
%  third. Repeated for each ratio.
%    A-C  leader fraction x(t), held-out 1:1, 1:2, 1:5 leader:follower
%    D-F  total cells N/N0 (log y-axis) for the same held-out co-cultures
%  Untreated = blue, +TGF-beta1 = orange.
%  Solid lines: prediction; shaded: 95% bootstrap prediction band (bootstrap
%  of the two training ratios); symbols: held-out data; dotted: fit to all ratios.
%
%  alpha and delta are read from calibrated_params.mat (Fig5 script);
%  if absent, the typed values below are used.
%  Output: ../figures/FigS9_leave_one_ratio_out.png / .pdf
%  Needs: Optimization + Statistics Toolboxes, MATLAB R2020a or newer.

clear; close all; clc; warning('off','all'); rng(1);

nBoot = 1000;                        % bootstrap samples per held-out ratio (100 for a test)
tData = (0:4)';
tGrid = linspace(0, 4, 101)';

%% ---------------- Style (same as Fig. 3) ----------------
baseFont = 'Helvetica';
LABEL_FS = 18;  TITLE_FS = 16;  TICK_FS = 11;  LEGEND_FS = 12;  PANEL_FS = 20;
FIT_LW = 2;  DATA_LW = 1.5;  MARKER_SIZE = 7;  BAND_ALPHA = 0.15;
condCol    = [0.00 0.45 0.74;        % untreated: blue
              0.85 0.33 0.10];       % +TGF-beta1: orange
condMk     = {'o', 's'};
condLabels = {'Untreated', '+TGF-\beta1'};
ratioLabels = {'1:1', '1:2', '1:5'}; % leader:follower

set(groot, 'DefaultFigureColor','w', 'DefaultAxesFontName',baseFont, ...
    'DefaultTextFontName',baseFont, 'DefaultAxesFontSize',TICK_FS, ...
    'DefaultAxesLabelFontSizeMultiplier',1, 'DefaultAxesTitleFontSizeMultiplier',1, ...
    'DefaultAxesLineWidth',1, 'DefaultAxesTickDir','out', 'DefaultAxesBox','off');

%% ---------------- Fixed growth rates ----------------
if isfile('calibrated_params.mat')
    cal = load('calibrated_params.mat');
    alphaFix = cal.alphaFix;  deltaFix = cal.deltaFix;
    fprintf('Loaded calibrated_params.mat\n');
else
    alphaFix = [0.3792 0.1212];      % untreated, +TGF-beta1
    deltaFix = [0.4470 0.4470];      % shared
end

%% ---------------- Co-culture data ----------------
% condition (1 untreated, 2 +TGF), ratio (1 = 1:1, 2 = 1:2, 3 = 1:5 L:F), leaders, total
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
    S(i).cond  = raw{i,1};  %#ok<SAGROW>
    S(i).ratio = raw{i,2};
    S(i).x     = L ./ N;
    S(i).n     = N ./ N(1);
    S(i).alpha = alphaFix(S(i).cond);
    S(i).delta = deltaFix(S(i).cond);
end
% residual scaling: pooled SDs over days 1-4 of ALL series (same as calibration)
allX = cell2mat(arrayfun(@(s) s.x(2:end), S(:), 'UniformOutput', false));
allN = cell2mat(arrayfun(@(s) s.n(2:end), S(:), 'UniformOutput', false));
W = [std(allX, 1), std(allN, 1)];

opts = optimoptions('lsqnonlin', 'Display','off', ...
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10, 'OptimalityTolerance',1e-10);

%% ---------------- Fit to all ratios (dotted reference) ----------------
bgAll = fit_shared(S, tData, W, opts);
fprintf('Fit to all ratios: beta = %.4f, gamma = %.4f\n', bgAll(1), bgAll(2));

%% ---------------- Leave-one-ratio-out ----------------
for r = 1:3
    train = find([S.ratio] ~= r);
    test  = find([S.ratio] == r);
    fprintf('Holding out %s L:F (%d bootstrap samples)...\n', ratioLabels{r}, nBoot);

    bgL = fit_shared(S(train), tData, W, opts);
    nTr = 2 * (numel(tData) - 1) * numel(train);
    c   = sqrt(nTr / (nTr - 2));
    PB  = nan(nBoot, 2);
    for b = 1:nBoot
        try
            PB(b,:) = fit_shared(boot_resample(S(train), bgL, tData, c), tData, W, opts);
        catch
        end
    end
    PB = PB(~any(isnan(PB), 2), :);

    LOO(r).bg      = bgL;                                   %#ok<SAGROW>
    LOO(r).bgCI    = prctile(PB, [2.5 97.5], 1);
    LOO(r).ssePred = sum(lf_resid(bgL,   S(test), tData, W).^2);
    LOO(r).sseAll  = sum(lf_resid(bgAll, S(test), tData, W).^2);
    for j = 1:numel(test)
        i = test(j);
        [xP, nP] = lf_sim(bgL,   S(i), tGrid);
        [xA, nA] = lf_sim(bgAll, S(i), tGrid);
        Xb = zeros(size(PB,1), numel(tGrid));  Nb = Xb;
        for b = 1:size(PB,1)
            [xx, nn] = lf_sim(PB(b,:), S(i), tGrid);
            Xb(b,:) = xx';  Nb(b,:) = nn';
        end
        LOO(r).ser(j) = struct('i',i, 'xP',xP, 'nP',nP, 'xA',xA, 'nA',nA, ...
            'xCI',prctile(Xb,[2.5 97.5],1), 'nCI',prctile(Nb,[2.5 97.5],1));
    end
end

fprintf('\n  Held out   beta    gamma    weighted SSE: prediction   fit to all\n');
for r = 1:3
    fprintf('  %-8s  %6.3f   %6.3f        %8.3f        %8.3f\n', ratioLabels{r}, ...
        LOO(r).bg(1), LOO(r).bg(2), LOO(r).ssePred, LOO(r).sseAll);
end

%% ---------------- Figure ----------------
fig = figure('Units','inches', 'Position',[1 1 14 8.5]);
tl  = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','loose');   % room for panel letters
letters = 'ABCDEF';
for row = 1:2                                   % 1 = x(t), 2 = N/N0
    for r = 1:3
        ax = nexttile(tl, (row-1)*3 + r);  hold(ax, 'on');
        set(ax, 'FontSize',TICK_FS, 'LineWidth',1, 'TickDir','out', 'Box','off', 'Layer','top');
        grid(ax, 'on');  ax.GridAlpha = 0.18;
        h = gobjects(1, 3);
        for j = 1:numel(LOO(r).ser)
            sr = LOO(r).ser(j);  i = sr.i;  cnd = S(i).cond;  col = condCol(cnd,:);
            if row == 1
                yP = sr.xP;  yC = sr.xCI;  yA = sr.xA;  yD = S(i).x;
            else
                yP = sr.nP;  yC = sr.nCI;  yA = sr.nA;  yD = S(i).n;
            end
            fill(ax, [tGrid; flipud(tGrid)], [yC(1,:)'; flipud(yC(2,:)')], col, ...
                 'FaceAlpha',BAND_ALPHA, 'EdgeColor','none', 'HandleVisibility','off');
            h(cnd) = plot(ax, tGrid, yP, '-', 'Color',col, 'LineWidth',FIT_LW);
            h(3)   = plot(ax, tGrid, yA, ':', 'Color',[0.2 0.2 0.2], 'LineWidth',DATA_LW);
            plot(ax, tData, yD, condMk{cnd}, 'Color',col, 'MarkerFaceColor',col, ...
                 'MarkerSize',MARKER_SIZE, 'LineWidth',DATA_LW, 'HandleVisibility','off');
        end
        if row == 1
            title(ax, sprintf('%s L:F Held Out', ratioLabels{r}), 'FontWeight','bold', 'FontSize',TITLE_FS);
            ylim(ax, [0 0.6]);
            if r == 1, ylabel(ax, 'Leader Fraction x(t)', 'FontWeight','bold', 'FontSize',LABEL_FS); end
        else
            set(ax, 'YScale','log');  ylim(ax, [0.9 9]);  set(ax, 'YTick',[1 2 4 8]);
            xlabel(ax, 'Time (days)', 'FontWeight','bold', 'FontSize',LABEL_FS);
            if r == 1, ylabel(ax, 'Total Cells (N/N_0)', 'FontWeight','bold', 'FontSize',LABEL_FS); end
        end
        if row == 1 && r == 1
            lg = legend(ax, h, {'Prediction, untreated', 'Prediction, +TGF-\beta1', 'Fit to all ratios'}, ...
                        'Location','northeast', 'Box','off');
            lg.FontSize = LEGEND_FS;  lg.ItemTokenSize = [20 10];
        end
        text(ax, -0.12, 1.04, letters((row-1)*3 + r), 'Units','normalized', 'FontSize',PANEL_FS, ...
             'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Clipping','off');
    end
end
drawnow;
outDir = fullfile('..', 'figures');  if ~isfolder(outDir), mkdir(outDir); end
exportgraphics(fig, fullfile(outDir, 'FigS9_leave_one_ratio_out.png'), 'Resolution',300);
exportgraphics(fig, fullfile(outDir, 'FigS9_leave_one_ratio_out.pdf'), 'ContentType','vector');
fprintf('\nSaved FigS9_leave_one_ratio_out.png / .pdf\n');

%% ================= Local functions =================
function dz = rhs_2D(~, z, p)
    beta = p(1);  gamma = p(2);  alpha = p(3);  delta = p(4);
    x = z(1);  N = z(2);
    A = alpha*x + beta*(1-x);
    B = gamma*x + delta*(1-x);
    dz = [(A - B)*x*(1-x); (A*x + B*(1-x))*N];
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
    bg = lsqnonlin(@(p) lf_resid(p, S, t, W), [0.3 0.3], [-2 -2], [3 3], opts);
end

function Sb = boot_resample(S, bg, t, c)
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
