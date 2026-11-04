%% ==========================================================
%  CROSS-PREDICTION: use 1:1 parameters only
%  Predict 2:1 and 5:1 without refitting
%% ==========================================================

warning('off','all');
rng(1);

tdata = (0:4)';
tspan = [tdata(1), tdata(end)];
tFine = linspace(tspan(1), tspan(2), 400)';

optsODE = odeset('RelTol',1e-9,'AbsTol',1e-12);

%% Data
D.FL21.name = '2:1 F:L';
D.FL21.y = [0.363; 0.31; 0.204705882; 0.171052632; 0.12];
D.FL21.N = [1.0E+06; 2.0E+06; 2.8E+06; 3.8E+06; 5.3E+06];

D.FL51.name = '5:1 F:L';
D.FL51.y = [0.178; 0.145945946; 0.086666667; 0.075; 0.054545455];
D.FL51.N = [9.9E+05; 1.9E+06; 3.0E+06; 4.4E+06; 7.3E+06];

predSets = {D.FL21, D.FL51};

%% Parameters from 1:1 only
beta_mean  = 0.126324094891;
gamma_mean = 0.728494058685;

beta_CI  = [-0.024640717587, 0.255249836867];
gamma_CI = [ 0.540035161686, 0.934653887606];

pMean = [beta_mean, gamma_mean];

%% Bootstrap-like uncertainty sampling inside CI rectangle
nDraw = 2000;
qLo = 2.5;
qHi = 97.5;

beta_draw  = beta_CI(1)  + (beta_CI(2)  - beta_CI(1))  * rand(nDraw,1);
gamma_draw = gamma_CI(1) + (gamma_CI(2) - gamma_CI(1)) * rand(nDraw,1);

%% Plot
figure('Units','inches','Position',[1 1 13 5]); clf;

sgtitle('Parameters from 1:1 predict 2:1 and 5:1', ...
    'FontSize', 12, 'FontWeight', 'bold');

cols = [
    0.85 0.33 0.10;
    0.00 0.45 0.74
];

%% Leader fraction prediction
subplot(1,2,1); hold on; box on;

for k = 1:2
    S = predSets{k};
    y0 = S.y(1);
    N0 = S.N(1);
    col = cols(k,:);

    [yMean, ~] = simulate_2D(pMean, tFine, y0, N0, optsODE);

    yMat = nan(nDraw, numel(tFine));

    for b = 1:nDraw
        [yb, ~] = simulate_2D([beta_draw(b), gamma_draw(b)], ...
                              tFine, y0, N0, optsODE);
        yMat(b,:) = yb(:).';
    end

    yMat = yMat(all(isfinite(yMat),2), :);

    yLo = prctile(yMat, qLo, 1);
    yHi = prctile(yMat, qHi, 1);

    fill([tFine; flipud(tFine)], ...
         [yLo(:); flipud(yHi(:))], ...
         col, 'FaceAlpha', 0.22, ...
         'EdgeColor','none', ...
         'HandleVisibility','off');

    plot(tFine, yMean, '-', ...
        'Color', col, ...
        'LineWidth', 2.6, ...
        'DisplayName', sprintf('%s prediction', S.name));

    plot(tdata, S.y, 'o', ...
        'Color', col, ...
        'MarkerFaceColor', col, ...
        'MarkerSize', 6, ...
        'LineWidth', 1.2, ...
        'HandleVisibility','off');
end

grid on;
set(gca,'Layer','top','FontSize',11);
xlabel('Time (days)','FontWeight','bold');
ylabel('Leader fraction y(t)','FontWeight','bold');
ylim([0 1]);
xlim(tspan);
legend('Location','best','Box','off');
title(sprintf('y(t): mean + %.1f--%.1f%% band', qLo, qHi), ...
    'FontSize', 11);

%% Total population prediction
subplot(1,2,2); hold on; box on;

for k = 1:2
    S = predSets{k};
    y0 = S.y(1);
    N0 = S.N(1);
    col = cols(k,:);

    [~, NMean] = simulate_2D(pMean, tFine, y0, N0, optsODE);

    NMat = nan(nDraw, numel(tFine));

    for b = 1:nDraw
        [~, Nb] = simulate_2D([beta_draw(b), gamma_draw(b)], ...
                              tFine, y0, N0, optsODE);
        NMat(b,:) = Nb(:).';
    end

    good = all(isfinite(NMat),2) & all(NMat > 0,2);
    NMat = NMat(good,:);

    NLo = prctile(NMat, qLo, 1);
    NHi = prctile(NMat, qHi, 1);

    fill([tFine; flipud(tFine)], ...
         [NLo(:); flipud(NHi(:))], ...
         col, 'FaceAlpha', 0.22, ...
         'EdgeColor','none', ...
         'HandleVisibility','off');

    plot(tFine, NMean, '-', ...
        'Color', col, ...
        'LineWidth', 2.6, ...
        'DisplayName', sprintf('%s prediction', S.name));

    plot(tdata, S.N, 'o', ...
        'Color', col, ...
        'MarkerFaceColor', col, ...
        'MarkerSize', 6, ...
        'LineWidth', 1.2, ...
        'HandleVisibility','off');
end

set(gca,'YScale','log');
grid on;
set(gca,'Layer','top','FontSize',11);
xlabel('Time (days)','FontWeight','bold');
ylabel('Total population N(t)','FontWeight','bold');
xlim(tspan);
legend('Location','best','Box','off');
title(sprintf('N(t): mean + %.1f--%.1f%% band', qLo, qHi), ...
    'FontSize', 11);

saveas(gcf, 'cross_prediction_from_1to1_only.png');

%% ==========================================================
% Local functions
%% ==========================================================

function [yOut, NOut] = simulate_2D(p, tvec, y0, N0, optsODE)

    try
        sol = ode23s(@(t,z) rhs_2D(t,z,p), ...
                     [tvec(1), tvec(end)], ...
                     [y0; N0], optsODE);

        Z = deval(sol, tvec);

        yOut = Z(1,:).';
        NOut = Z(2,:).';

        if any(~isfinite(yOut)) || any(~isfinite(NOut)) || any(NOut <= 0)
            yOut = nan(size(tvec));
            NOut = nan(size(tvec));
        end

    catch
        yOut = nan(size(tvec));
        NOut = nan(size(tvec));
    end
end

function dzdt = rhs_2D(~, z, p)

    beta  = p(1);
    gamma = p(2);

    y = z(1);
    N = z(2);

    A = 0.3798*y + beta*(1-y);
    B = gamma*y  + 0.4467*(1-y);

    dydt = (A-B)*y*(1-y);
    dNdt = (A*y + B*(1-y))*N;

    dzdt = [dydt; dNdt];
end