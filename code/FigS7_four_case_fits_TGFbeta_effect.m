%% Supplementary Fig. S7: fits of the four cases for the effect of TGF-beta1 (label fig:hypfits)
%
%  The co-culture model is fitted jointly to all six co-cultures (1:1, 1:2,
%  1:5 leader:follower; untreated and +TGF-beta1) under four cases:
%     H0      beta and gamma shared by both conditions      (black)
%     Hbeta   only beta differs between conditions          (blue)
%     Hgamma  only gamma differs between conditions         (orange)
%     Hbg     both differ                                   (purple)
%  alpha and delta are fixed from the monoculture fits.
%    A-B  leader fraction x(t), untreated (A) and +TGF-beta1 (B)
%    C-D  total cells N/N0 (log y-axis), untreated (C) and +TGF-beta1 (D)
%  Color = case, line style + marker = seeding ratio.
%  Because the four fits nearly overlap, they are drawn as nested lines:
%  H0 widest underneath, the other cases progressively thinner on top.
%
%  Reads alpha, delta from calibrated_params.mat (Fig5 script) if present.
%  Output: ../figures/FigS7_four_case_fits.png / .pdf. Runs in seconds (no bootstrap).

clear; close all; clc; warning('off','all');

tData = (0:4)';
tGrid = linspace(0, 4, 101)';

%% ---------------- Style (same as Fig. 5) ----------------
baseFont = 'Helvetica';
LABEL_FS = 18;  TITLE_FS = 16;  TICK_FS = 11;  LEGEND_FS = 12;  PANEL_FS = 20;
FIT_LW = 2;  DATA_LW = 1.5;  MARKER_SIZE = 7;
caseCol  = [0.00 0.00 0.00;          % H0      black
            0.00 0.45 0.74;          % Hbeta   blue
            0.85 0.33 0.10;          % Hgamma  orange
            0.49 0.18 0.56];         % Hbg     purple
caseLW   = [5.0 3.4 2.2 1.2];        % nested widths so overlapping fits stay visible
ratioLS  = {'-', '--', ':'};         % 1:1 solid, 1:2 dashed, 1:5 dotted
ratioMk  = {'o', 's', '^'};
dataCol  = [0.35 0.35 0.35];         % data markers (gray), shape = ratio
ratioLabels = {'1:1 L:F', '1:2 L:F', '1:5 L:F'};
condLabels  = {'Untreated', '+TGF-\beta1'};
hypNames = {'H_0: neither changes', 'H_\beta: only \beta changes', ...
            'H_\gamma: only \gamma changes', 'H_{\beta\gamma}: both change'};
hypK     = [2 3 3 4];

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
allX = cell2mat(arrayfun(@(s) s.x(2:end), S(:), 'UniformOutput', false));
allN = cell2mat(arrayfun(@(s) s.n(2:end), S(:), 'UniformOutput', false));
W = [std(allX, 1), std(allN, 1)];
nRes = 2 * (numel(tData) - 1) * nSer;          % 48

opts = optimoptions('lsqnonlin', 'Display','off', ...
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10, 'OptimalityTolerance',1e-10);

%% ---------------- Fit the four cases ----------------
q = zeros(4, 4);  sse = zeros(1, 4);  aic = zeros(1, 4);
for h = 1:4
    k   = hypK(h);
    fun = @(p) case_resid(expand_params(p, h), S, tData, W);
    p   = lsqnonlin(fun, 0.3*ones(1,k), -2*ones(1,k), 3*ones(1,k), opts);
    q(h,:) = expand_params(p, h);              % [beta_U gamma_U beta_T gamma_T]
    sse(h) = sum(fun(p).^2);
    kk = k + 1;
    aic(h) = nRes*log(sse(h)/nRes) + 2*kk + 2*kk*(kk+1)/(nRes - kk - 1);
end
fprintf('\n  Case      beta_U  gamma_U  beta_T  gamma_T   SSE     dAICc\n');
for h = 1:4
    fprintf('  %-8s  %6.3f  %6.3f   %6.3f  %6.3f   %6.3f   %5.2f\n', ...
        strtok(strrep(strrep(hypNames{h},'\',''),'_',''),':'), q(h,:), sse(h), aic(h)-min(aic));
end

%% ---------------- Figure ----------------
fig = figure('Units','inches', 'Position',[1 1 11 9.5]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','loose');
letters = 'ABCD';
for row = 1:2                                  % 1 = x(t), 2 = N/N0
    for cnd = 1:2
        tileNo = (row-1)*2 + cnd;
        ax = nexttile(tl, tileNo);  hold(ax, 'on');
        set(ax, 'FontSize',TICK_FS, 'LineWidth',1, 'TickDir','out', 'Box','off', 'Layer','top');
        grid(ax, 'on');  ax.GridAlpha = 0.18;
        cols = 2*cnd - 1 : 2*cnd;
        for i = find([S.cond] == cnd)
            rr = S(i).ratio;
            for h = 1:4                                % H0 first (widest), then thinner on top
                [xx, nn] = lf_sim(q(h, cols), S(i), tGrid);
                if row == 1, yy = xx; else, yy = nn; end
                plot(ax, tGrid, yy, ratioLS{rr}, 'Color',caseCol(h,:), 'LineWidth',caseLW(h), ...
                     'HandleVisibility','off');
            end
            if row == 1, yD = S(i).x; else, yD = S(i).n; end
            plot(ax, tData, yD, ratioMk{rr}, 'Color',dataCol, 'MarkerFaceColor',dataCol, ...
                 'MarkerEdgeColor','w', 'MarkerSize',MARKER_SIZE+1, 'LineWidth',0.75, ...
                 'HandleVisibility','off');
        end
        if row == 1
            title(ax, ['Co-culture, ' condLabels{cnd}], 'FontWeight','bold', 'FontSize',TITLE_FS);
            ylabel(ax, 'Leader Fraction x(t)', 'FontWeight','bold', 'FontSize',LABEL_FS);
            ylim(ax, [0 0.6]);
        else
            set(ax, 'YScale','log');  ylim(ax, [0.9 9]);  set(ax, 'YTick',[1 2 4 8]);
            ylabel(ax, 'Total Cells (N/N_0)', 'FontWeight','bold', 'FontSize',LABEL_FS);
            xlabel(ax, 'Time (days)', 'FontWeight','bold', 'FontSize',LABEL_FS);
        end
        % legends: color = case (panel A), line style + marker = ratio (panel B)
        if tileNo == 1
            hh = gobjects(1, 4);
            for h = 1:4
                hh(h) = plot(ax, NaN, NaN, '-', 'Color',caseCol(h,:), 'LineWidth',max(caseLW(h),2.5));
            end
            lg = legend(ax, hh, hypNames, 'Location','northeast', 'Box','off', 'Interpreter','tex');
            lg.FontSize = LEGEND_FS;  lg.ItemTokenSize = [28 10];
        elseif tileNo == 2
            hr = gobjects(1, 3);
            for rr = 1:3
                hr(rr) = plot(ax, NaN, NaN, [ratioLS{rr} ratioMk{rr}], 'Color',[0.2 0.2 0.2], ...
                    'MarkerFaceColor',dataCol, 'MarkerEdgeColor','w', 'MarkerSize',MARKER_SIZE+1, ...
                    'LineWidth',FIT_LW);
            end
            lg = legend(ax, hr, ratioLabels, 'Location','northeast', 'Box','off');
            lg.FontSize = LEGEND_FS;  lg.ItemTokenSize = [32 10];
        end
        text(ax, -0.12, 1.04, letters(tileNo), 'Units','normalized', 'FontSize',PANEL_FS, ...
             'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Clipping','off');
    end
end
drawnow;
outDir = fullfile('..', 'figures');  if ~isfolder(outDir), mkdir(outDir); end
exportgraphics(fig, fullfile(outDir, 'FigS7_four_case_fits.png'), 'Resolution',300);
exportgraphics(fig, fullfile(outDir, 'FigS7_four_case_fits.pdf'), 'ContentType','vector');
fprintf('\nSaved FigS7_four_case_fits.png / .pdf\n');

%% ================= Local functions =================
function q = expand_params(p, h)
% [beta_U gamma_U beta_T gamma_T] for each case
    switch h
        case 1, q = [p(1) p(2) p(1) p(2)];      % H0
        case 2, q = [p(1) p(2) p(3) p(2)];      % only beta changes
        case 3, q = [p(1) p(2) p(1) p(3)];      % only gamma changes
        otherwise, q = p(1:4);                  % both change
    end
end

function r = case_resid(q, S, t, W)
    r = [];
    for i = 1:numel(S)
        bg = q(2*S(i).cond - 1 : 2*S(i).cond);
        [x, n] = lf_sim(bg, S(i), t);
        r = [r; (x(2:end) - S(i).x(2:end)) / W(1);
                (n(2:end) - S(i).n(2:end)) / W(2)]; %#ok<AGROW>
    end
end

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
