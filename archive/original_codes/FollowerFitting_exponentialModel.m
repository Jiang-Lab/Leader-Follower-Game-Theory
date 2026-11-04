%Final Follower Code

% Load tumorall; clc; format long; warning('off','all'); rng(1);  growth data (time series in 6 replicates)
tdata = [0; 1; 2; 3; 4];  % Time points
Tdata = [1,	1,		1,	1,	1;
1.257692308,	1.594160584,			1.505535055,	1.52398524,	1.8E+00;
1.729807692	,1.052554745,		1.236162362	,1.110701107	, 2.8E+00
3.576923077,	2.832116788	,		3.837638376	,3.837638376,	4.2E+00;
6.6	, 4.986131387,	6.715867159,	6.3099631,	6.7E+00];

% SD across replicates at each time (for error bars you provided)
st = [0, 0.171619374,0.633706604,0.44980746,0.650221324];

% Bounds + initial guess for parameters [growth_rate, initial_size]
lb = [0.0, min(Tdata(:))];
ub = [1.0, max(Tdata(:))];
x0 = [0.5, Tdata(1,1)];

% Optimizer settings for lsqcurvefit
options = optimoptions('lsqcurvefit', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12);

% ---------- Bootstrap (block-resample columns + percentile CI) ----------
nBoot = 10;                                  % use a few thousand for stable quantiles
[~, nRep] = size(Tdata);
paramEsts = zeros(nBoot, numel(x0));

for b = 1:nBoot
    % 1) resample replicates as blocks (columns)
    idx = randsample(nRep, nRep, true);
    Tb  = Tdata(:, idx);

    % 2) target we fit: mean across resampled replicates
    syntheticData = mean(Tb, 2);

    % 3) fit model to this bootstrap sample
    try
        xBoot = lsqcurvefit(@(x,t) cost_function(x,t), x0, tdata, syntheticData, lb, ub, options);
        paramEsts(b,:) = xBoot(:).';
    catch
        paramEsts(b,:) = NaN;
    end
end

% Remove failed fits (if any)
paramEsts = paramEsts(~any(isnan(paramEsts),2), :);

% Parameter summary and 95% percentile CI (NO /sqrt(nBoot))
paramMean = mean(paramEsts,1);
paramCI   = squeeze(prctile(paramEsts,[2.5 97.5],1));   % 2 x numParams

% Prediction band from bootstrap parameter draws
tgrid = linspace(min(tdata), max(tdata), 200);
yBoot = zeros(size(paramEsts,1), numel(tgrid));
for b = 1:size(paramEsts,1)
    yBoot(b,:) = cost_function(paramEsts(b,:), tgrid);
end
predMean = mean(yBoot,1);
predCI   = squeeze(prctile(yBoot,[2.5 97.5],1));        % 2 x length(tgrid)

% Display parameter estimates and percentile CIs
disp('Mean of Bootstrapped Parameters:');
disp(paramMean);
disp('95% Percentile CIs (rows: [low; high])');
disp(paramCI);

%% ---------- PNAS-like style (session defaults; no extra functions) ----------
baseFont = 'Helvetica';            % or 'Arial'
PNAS_LABEL_FONTSIZE  = 18;
PNAS_TICK_FONTSIZE   = 11;
PNAS_LEGEND_FONTSIZE = 12;

% Color order (magenta, green, blue, red-orange, orange)
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

%% ---------- Plot (uniform fonts, legend, colors, bold labels) ----------
figure('Units','inches','Position',[1 1 4.5 3.5]); hold on;

% Experimental data with error bars
errorbar(tdata, mean(Tdata,2), st(:), 'k*', ...
    'LineWidth',1.5, 'MarkerSize',6, 'CapSize',6, ...
    'DisplayName','Mean Data');

% Shaded 95% prediction band
fill([tgrid, fliplr(tgrid)], [predCI(1,:), fliplr(predCI(2,:))], ...
     [0 0 0], 'FaceAlpha',0.08, 'EdgeColor','none', 'DisplayName','95% band');

% Mean fit (from bootstrap mean params OR from predMean)
plot(tgrid, predMean, 'k-', 'LineWidth',2, 'DisplayName','Mean Fit');

% Axis labels (bold, uniform size)
xlabel('Time (days)',   'FontWeight','bold', 'FontSize',PNAS_LABEL_FONTSIZE);
ylabel('Follower Cells','FontWeight','bold', 'FontSize',PNAS_LABEL_FONTSIZE);

% Axes/legend aesthetics + panel letter outside upper-left
ax = gca; grid(ax,'on'); ax.Layer = 'top';
lg = legend('Location','eastoutside','Box','off'); 
lg.FontSize = PNAS_LEGEND_FONTSIZE; lg.ItemTokenSize = [16 10];

text(ax, -0.12, 1.04, 'A', 'Units','normalized', ...
    'FontSize',20, 'FontWeight', ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
    'Clipping','off');

% Save figure
saveas(gcf, 'model_fit_with_CI_PNASstyle.png');
% exportgraphics(gcf, 'model_fit_with_CI_PNASstyle.pdf', 'ContentType','vector');

% ----------------- Model (exponential growth) -----------------
function dydt = g(t, y, x)
    % x(1) = growth rate, x(2) = initial size
    dydt = (x(1)) * y;
end

function y = cost_function(x, tdata)
    % Solve y' = x1*y with y(0) = x2 at the first t in tdata
    [~, Y] = ode23s(@(t,y) g(t,y,x), tdata, x(2));
    y = Y(:,1);
end
