%% ================================================================
%  Sensitivity analysis for 2D system (y, N) with TWO datasets
%  LOCAL OAT (±p%) + GLOBAL Monte Carlo (±p% jointly)
%  Model:
%     A = alpha*y + beta*(1-y)
%     B = gamma*y + delta*(1-y)
%     y' = (A - B) * y * (1 - y)
%     fbar = A*y + B*(1-y)
%     N' = fbar * N
%
%  Plots include:
%   - Base fits (y and N)
%   - Local OAT sensitivity bars for y and N separately
%   - Timewise local deviation curves |Δy(t)| and |ΔN(t)|
%   - Global game robustness + distribution
%   - Global SRC + PRCC for multiple outputs including RMSEy and RMSEn
%   - beta vs delta boundary scatter
%
%  Author: you + ChatGPT
%  ================================================================

 format long; warning('off','all');

%% ---------------------- Data & Base Setup -----------------------
tdata = [0; 1; 2; 3; 4];

yobs = [0.521; 0.447058824; 0.3375; 0.294117647; 0.216923077];
Nobs = [8.7E+05; 1.9E+06; 2.5E+06; 2.8E+06; 3.5E+06];

% Initial conditions (use first observation)
y0 = yobs(1);
N0 = Nobs(1);
z0 = [y0; N0];

param_names = {'\alpha','\beta','\gamma','\delta'};
base_params = [0.3798, 0.1263, 0.7285, 0.4467];  % [alpha, beta, gamma, delta]

perturb_percent_local = 0.10;   % ±10% local OAT
range_global = 0.25;            % ±25% global joint range
Nmc = 5000;                     % Monte Carlo samples

% weights for joint metric (used for RMSEjoint only)
wY = 1.0;
wN = 1.0;

%% ---------------------- Solver Options --------------------------
opts = odeset('NonNegative',1:2,'RelTol',1e-8,'AbsTol',1e-10);

%% ---------------------- Base Simulation -------------------------
[~, Zbase] = ode23s(@(t,z) g_rhs_2d(t,z,base_params), tdata, z0, opts);
ybase = Zbase(:,1);
Nbase = Zbase(:,2);

base_game = classify_game(base_params);
fprintf('Base game: %s\n', base_game);

rmse_y_base = rmse(ybase, yobs);
rmse_N_base = rmse(Nbase, Nobs);
fprintf('Base RMSE: y=%.6g, N=%.6g\n', rmse_y_base, rmse_N_base);

% Plot base trajectories vs data
figure('Name','Base fit: y(t)','Color','w');
plot(tdata, yobs, 'o', 'LineWidth',1.5); hold on; grid on; box on;
plot(tdata, ybase,'-','LineWidth',2);
xlabel('Time'); ylabel('y');
legend('y_{obs}','y_{model}','Location','best');
title('Base trajectory: y(t)');
hold off;

figure('Name','Base fit: N(t)','Color','w');
plot(tdata, Nobs, 'o', 'LineWidth',1.5); hold on; grid on; box on;
plot(tdata, Nbase,'-','LineWidth',2);
xlabel('Time'); ylabel('N');
legend('N_{obs}','N_{model}','Location','best');
title('Base trajectory: N(t)');
hold off;

%% ================================================================
%  PART A. LOCAL (OAT) SENSITIVITY  (y and N separately)
%  SIy = || y_pert - y_base || / || y_base ||
%  SIN = || N_pert - N_base || / || N_base ||
%  Also timewise deviations: |Δy(t)| and |ΔN(t)|
%  ================================================================
nParams = numel(base_params);
sens_local_y = zeros(nParams,2);  % [down, up]
sens_local_N = zeros(nParams,2);  % [down, up]

Dy_down = zeros(numel(tdata), nParams);
Dy_up   = zeros(numel(tdata), nParams);
DN_down = zeros(numel(tdata), nParams);
DN_up   = zeros(numel(tdata), nParams);

% Optional: RMSE-based local sensitivity (relative change in RMSE)
scaleN = mean(Nobs);  % normalize N residuals for joint RMSE
rmse_joint_base = sqrt( wY*mean((ybase - yobs).^2) + wN*mean(((Nbase - Nobs)/scaleN).^2) );
sens_local_rmse_joint = zeros(nParams,2);

games_local_down = strings(nParams,1);
games_local_up   = strings(nParams,1);

for i = 1:nParams
    p_down = base_params; p_up = base_params;
    p_down(i) = base_params(i) * (1 - perturb_percent_local);
    p_up(i)   = base_params(i) * (1 + perturb_percent_local);

    [~, Zd] = ode23s(@(t,z) g_rhs_2d(t,z,p_down), tdata, z0, opts);
    [~, Zu] = ode23s(@(t,z) g_rhs_2d(t,z,p_up),   tdata, z0, opts);

    yd = Zd(:,1); Nd = Zd(:,2);
    yu = Zu(:,1); Nu = Zu(:,2);

    % Separate trajectory sensitivity indices
    sens_local_y(i,1) = norm(yd - ybase,2) / max(norm(ybase,2),eps);
    sens_local_y(i,2) = norm(yu - ybase,2) / max(norm(ybase,2),eps);

    sens_local_N(i,1) = norm(Nd - Nbase,2) / max(norm(Nbase,2),eps);
    sens_local_N(i,2) = norm(Nu - Nbase,2) / max(norm(Nbase,2),eps);

    % Timewise deviations
    Dy_down(:,i) = abs(yd - ybase);
    Dy_up(:,i)   = abs(yu - ybase);

    DN_down(:,i) = abs(Nd - Nbase);
    DN_up(:,i)   = abs(Nu - Nbase);

    % Local RMSE-joint sensitivity (relative)
    rmse_joint_d = sqrt( wY*mean((yd - yobs).^2) + wN*mean(((Nd - Nobs)/scaleN).^2) );
    rmse_joint_u = sqrt( wY*mean((yu - yobs).^2) + wN*mean(((Nu - Nobs)/scaleN).^2) );
    sens_local_rmse_joint(i,1) = (rmse_joint_d - rmse_joint_base) / max(rmse_joint_base,eps);
    sens_local_rmse_joint(i,2) = (rmse_joint_u - rmse_joint_base) / max(rmse_joint_base,eps);

    % Game flips
    games_local_down(i) = classify_game(p_down);
    games_local_up(i)   = classify_game(p_up);
end

% -------- BAR: local sensitivity for y --------
figure('Name','Local OAT sensitivity: y only','Color','w');
bar(sens_local_y,'grouped'); grid on;
set(gca,'XTickLabel',param_names,'FontSize',12);
ylabel('Sensitivity Index');
legend(sprintf('-%d%%',round(100*perturb_percent_local)), sprintf('+%d%%',round(100*perturb_percent_local)), ...
       'Location','northwest');
annotation('textbox',[0.02 0.92 0.06 0.08],'String','B', ...
    'FontWeight','bold','FontSize',18,'LineStyle','none');
title(sprintf('Local OAT sensitivity of y(t) (±%d%%)', round(100*perturb_percent_local)));

% -------- BAR: local sensitivity for N --------
figure('Name','Local OAT sensitivity: N only','Color','w');
bar(sens_local_N,'grouped'); grid on;
set(gca,'XTickLabel',param_names,'FontSize',12);
ylabel('Sensitivity Index');
ylim([0 0.3]);
annotation('textbox',[0.02 0.92 0.06 0.08],'String','B', ...
    'FontWeight','bold','FontSize',18,'LineStyle','none');
legend(sprintf('-%d%%',round(100*perturb_percent_local)), sprintf('+%d%%',round(100*perturb_percent_local)), ...
       'Location','northwest');
title(sprintf('Local OAT sensitivity of N(t) (±%d%%)', round(100*perturb_percent_local)));

% -------- BAR: local sensitivity for joint RMSE (optional) --------
figure('Name','Local OAT sensitivity: joint RMSE change','Color','w');
bar(sens_local_rmse_joint,'grouped'); grid on;
set(gca,'XTickLabel',param_names,'FontSize',12);
ylabel('Relative change in RMSE_{joint}');
legend(sprintf('-%d%%',round(100*perturb_percent_local)), sprintf('+%d%%',round(100*perturb_percent_local)), ...
       'Location','northwest');
title(sprintf('Local effect on RMSE_{joint} (±%d%%)', round(100*perturb_percent_local)));

% -------- TIMEWISE: |Δy(t)| for + perturbation --------
figure('Name','Timewise local deviation: |Δy(t)| (+)','Color','w'); hold on; grid on; box on;
for i=1:nParams
    plot(tdata, Dy_up(:,i), '-o', 'LineWidth',1.5);
end
xlabel('Time'); ylabel('|y_{+} - y_{base}|');
legend(param_names,'Location','best');
title(sprintf('Timewise |Δy(t)| for +%d%% perturbation', round(100*perturb_percent_local)));
hold off;

% -------- TIMEWISE: |ΔN(t)| for + perturbation --------
figure('Name','Timewise local deviation: |ΔN(t)| (+)','Color','w'); hold on; grid on; box on;
for i=1:nParams
    plot(tdata, DN_up(:,i), '-o', 'LineWidth',1.5);
end
xlabel('Time'); ylabel('|N_{+} - N_{base}|');
legend(param_names,'Location','best');
title(sprintf('Timewise |ΔN(t)| for +%d%% perturbation', round(100*perturb_percent_local)));
hold off;

% --- Print local game flips
fprintf('\nLocal OAT game outcomes (±%.0f%%):\n', 100*perturb_percent_local);
for i=1:nParams
    fprintf('  %s: - -> %s | + -> %s\n', param_names{i}, games_local_down(i), games_local_up(i));
end

%% ================================================================
%  PART B. GLOBAL SENSITIVITY (Monte Carlo, ±range_global jointly)
%  Collect outputs for sensitivity:
%    y_end, N_end, RMSEy, RMSEn, RMSEjoint
%  Compute:
%    SRC (linear), PRCC (rank-based)
%  Game robustness & distribution
%  ================================================================
rng(1);
mults = (1 - range_global) + (2*range_global)*rand(Nmc,4);
P = base_params .* mults;

Yend  = zeros(Nmc,1);
Nend  = zeros(Nmc,1);
Ynrm  = zeros(Nmc,1);
Nnrm  = zeros(Nmc,1);
RMSEy = zeros(Nmc,1);
RMSEn = zeros(Nmc,1);
RMSEj = zeros(Nmc,1);
codes = zeros(Nmc,1);

for k = 1:Nmc
    p = P(k,:);
    [~, Z] = ode23s(@(t,z) g_rhs_2d(t,z,p), tdata, z0, opts);
    y = Z(:,1);
    N = Z(:,2);

    Yend(k) = y(end);
    Nend(k) = N(end);
    Ynrm(k) = norm(y,2);
    Nnrm(k) = norm(N,2);

    RMSEy(k) = rmse(y, yobs);
    RMSEn(k) = rmse(N, Nobs);
    RMSEj(k) = sqrt( wY*mean((y - yobs).^2) + wN*mean(((N - Nobs)/scaleN).^2) );

    codes(k) = game_code(p);
end

% --- Game robustness
base_code = game_code(base_params);
base_game = classify_game(base_params);
stay_frac = mean(codes == base_code);
flip_frac = 1 - stay_frac;

flip_to_dist = accumarray(codes(codes~=base_code), 1, [4,1]);
if any(codes~=base_code)
    flip_to_dist = flip_to_dist / sum(codes~=base_code);
else
    flip_to_dist = zeros(4,1);
end

fprintf('\nGLOBAL (±%.0f%% jointly) game robustness:\n', 100*range_global);
fprintf('  Stay in base game (%s): %.1f%%\n', base_game, 100*stay_frac);
fprintf('  Flip overall: %.1f%%\n', 100*flip_frac);
fprintf('  Flip to -> SnowDrift: %.1f%% | Harmony: %.1f%% | Coordination: %.1f%% | PD: %.1f%%\n', ...
    100*flip_to_dist(1), 100*flip_to_dist(2), 100*flip_to_dist(3), 100*flip_to_dist(4));

% --- Game distribution plot
figure('Name','Game type under global perturbations','Color','w');
cats = categorical({'Snow Drift','Harmony','Coordination','Prisoner''s Dilemma'});
counts = [sum(codes==1), sum(codes==2), sum(codes==3), sum(codes==4)];
bar(cats, counts); grid on;
ylabel('Count');
title(sprintf('Game types from global ±%.0f%% sampling', 100*range_global));

%% ---------------------- Global: SRC -----------------------------
outputs = {'Yend','Nend','RMSEjoint','RMSEy','RMSEn'};
Ymat = [Yend, Nend, RMSEj, RMSEy, RMSEn];

Zp = (P - mean(P,1)) ./ std(P,[],1);  % standardized params
SRC = zeros(4, size(Ymat,2));

for j = 1:size(Ymat,2)
    yj = Ymat(:,j);
    yj_std = (yj - mean(yj)) / std(yj);
    B = (Zp' * Zp) \ (Zp' * yj_std);
    SRC(:,j) = B;
end

figure('Name','SRC (Global sensitivity)','Color','w');
bar(SRC'); grid on;
set(gca,'XTickLabel',outputs,'FontSize',11);
legend(param_names,'Location','best');
ylabel('SRC');
title('Global sensitivity via SRC (higher |SRC| => stronger influence)');

%% ---------------------- Global: PRCC ----------------------------
PRCC = zeros(4, size(Ymat,2));
for j = 1:size(Ymat,2)
    PRCC(:,j) = partial_rank_corr(Ymat(:,j), P);
end

figure('Name','PRCC (Global sensitivity)','Color','w');
bar(PRCC'); grid on;
set(gca,'XTickLabel',outputs,'FontSize',11);
legend(param_names,'Location','best');
ylabel('PRCC');
title('Global sensitivity via PRCC (monotone/nonlinear)');

%% ---------------------- beta vs delta boundary scatter -----------
figure('Name','\beta vs \delta with boundary','Color','w'); hold on; grid on; box on;
b = P(:,2); d = P(:,4);
scatter(b(d>=b), d(d>=b), 8, 'filled');
scatter(b(d<b),   d(d<b),  8, 'filled');
plot([min(b) max(b)], [min(b) max(b)], 'k--', 'LineWidth',1);
xlabel('\beta'); ylabel('\delta');
legend('\delta \ge \beta','\delta < \beta','\beta=\delta','Location','best');
title('\beta vs \delta under global sampling');
hold off;

%% ---------------------- boundary margins -------------------------
beta  = base_params(2);
delta = base_params(4);
p_up_beta    = (delta/beta - 1) * 100;
p_down_delta = (1 - beta/delta) * 100;
p_sym        = ((delta - beta)/(delta + beta)) * 100;

fprintf('\nBoundary margins from original values (beta=%.6g, delta=%.6g):\n', beta, delta);
fprintf('  +%% needed on beta only to hit beta=delta:  %.6f%%\n', p_up_beta);
fprintf('  -%% needed on delta only to hit beta=delta: %.6f%%\n', p_down_delta);
fprintf('  Symmetric (beta up & delta down) to boundary: %.6f%%\n', p_sym);

%% =========================== FUNCTIONS ====================================

function dzdt = g_rhs_2d(~, z, x)
    % z = [y; N], x = [alpha, beta, gamma, delta]
    y = z(1); N = z(2);

    a = x(1); b = x(2); g = x(3); d = x(4);

    A = a*y + b*(1-y);
    B = g*y + d*(1-y);

    dydt = (A - B) * y * (1 - y);

    fbar = A*y + B*(1-y);
    dNdt = fbar * N;

    dzdt = [dydt; dNdt];
end

function lab = classify_game(p)
    a=p(1); b=p(2); g=p(3); d=p(4);
    if  (b>d && g>a), lab = "Snow Drift";
    elseif (b>d && g<a), lab = "Harmony";
    elseif (b<d && g<a), lab = "Coordination";
    else, lab = "Prisoner''s Dilemma";
    end
end

function code = game_code(p)
    % 1=SD, 2=Harmony, 3=Coordination, 4=PD
    a=p(1); b=p(2); g=p(3); d=p(4);
    if  (b>d && g>a), code = 1;
    elseif (b>d && g<a), code = 2;
    elseif (b<d && g<a), code = 3;
    else, code = 4;
    end
end

function val = rmse(x, y)
    x = x(:); y = y(:);
    val = sqrt(mean((x - y).^2));
end

function prcc_params = partial_rank_corr(Y, X)
    % PRCC between Y and each column of X, controlling for other X columns.
    % Y: N x 1, X: N x P
    P = size(X,2);
    rX = tiedrank(X);
    rY = tiedrank(Y(:));

    prcc_params = zeros(P,1);
    for j = 1:P
        others = setdiff(1:P, j);

        Bj   = (rX(:,others)' * rX(:,others)) \ (rX(:,others)' * rX(:,j));
        resX = rX(:,j) - rX(:,others) * Bj;

        By   = (rX(:,others)' * rX(:,others)) \ (rX(:,others)' * rY);
        resY = rY - rX(:,others) * By;

        prcc_params(j) = corr(resX, resY, 'type','Pearson');
    end
end