%% Table 3 and Table S2: fitted parameters, interior point x*, game and outcome for both calibrations
%
%  Including the first day (main fit): data normalized to day 0, days 0-4 fitted.
%  Excluding the first day: data normalized to day 1, days 1-4 fitted.
%
%  For each calibration:
%    monocultures : y = y0*exp(r t), replicate bootstrap -> alpha (untreated, +TGF-beta1),
%                   delta (one shared value, untreated estimate)
%    co-cultures  : the four cases H0, Hbeta, Hgamma, Hbg fitted jointly to all six
%                   co-cultures (leader fraction x and N/N0), alpha and delta fixed;
%                   AICc, residual-bootstrap 95% CIs, interior point
%                   x* = (beta - delta) / ((beta - delta) + (gamma - alpha)),
%                   game class, share of bootstrap fits in the same class, and the stable
%                   outcome for the observed initial leader fractions.
%
%  The day-0 part reproduces the Fig5 script exactly (same rng(1) sequence for the
%  monocultures and the H0 bootstrap), so Table 3 and Fig. 5 report the same numbers.
%  The day-1 part starts from rng(2), as in the Fig5 script (panels G-H).
%
%  Output: table rows printed to the Command Window, and
%          calibrated_params_both.mat (all estimates, both calibrations).
%  Runtime: about 30 min with nBoot = 1000 (eight bootstraps). nBoot = 100 for a quick test.
%  Needs: Optimization + Statistics Toolboxes, MATLAB R2020a or newer.

clear; close all; clc; warning('off','all');

nBoot = 1000;

%% ---------------- Data ----------------
% Monocultures normalized to day 0; rows = days 0-4, columns = 5 replicates
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
% Co-cultures: condition (1 untreated, 2 +TGF), ratio (1 = 1:1, 2 = 1:2, 3 = 1:5 L:F), leaders, total
raw = { ...
 1, 1, [4.5e5 7.6e5 9.0e5 1.0e6 9.4e5], [8.7e5 1.7e6 2.7e6 3.4e6 4.3e6];
 1, 2, [3.6e5 6.2e5 5.8e5 6.5e5 6.4e5], [1.0e6 2.0e6 2.8e6 3.8e6 5.3e6];
 1, 3, [1.8e5 2.7e5 2.6e5 3.3e5 4.0e5], [9.9e5 1.9e6 3.0e6 4.4e6 7.3e6];
 2, 1, [4.5e5 8.1e5 7.8e5 6.6e5 5.3e5], [8.7e5 1.9e6 2.5e6 2.8e6 3.5e6];
 2, 2, [3.6e5 5.4e5 5.5e5 5.5e5 4.3e5], [1.0e6 2.0e6 3.0e6 4.0e6 5.7e6];
 2, 3, [1.8e5 2.6e5 2.4e5 2.7e5 2.3e5], [9.9e5 2.0e6 3.0e6 4.4e6 6.7e6]};

caseNames = {'H0', 'Hbeta', 'Hgamma', 'Hbg'};
caseK     = [2 3 3 4];
condTxt   = {'untreated', '+TGF-beta1'};
optsMono = optimoptions('lsqcurvefit', 'Display','off', ...
    'FunctionTolerance',1e-12, 'StepTolerance',1e-12, 'OptimalityTolerance',1e-12);
opts = optimoptions('lsqnonlin', 'Display','off', ...
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10, 'OptimalityTolerance',1e-10);

startDay = [0 1];
for cal = 1:2
    s0 = startDay(cal);
    rng(cal);                                   % rng(1) for day 0 (as Fig5), rng(2) for day 1
    fitName = {'Including the first day', 'Excluding the first day'};
    fprintf('\n================ %s (days %d-4, normalized to day %d) ================\n', fitName{cal}, s0, s0);

    %% ---- monocultures ----
    tM = (0:4-s0)';
    for ct = 1:2
        for cnd = 1:2
            T = mono{ct,cnd}(1+s0:end, :);
            T = T ./ T(1, :);                   % normalize each replicate to the start day
            M(ct,cnd) = mono_fit(T, tM, nBoot, optsMono); %#ok<SAGROW>
        end
    end
    alphaFix = [M(1,1).rMean, M(1,2).rMean];
    deltaFix = [M(2,1).rMean, M(2,1).rMean];   % one shared follower rate (untreated)
    fprintf('Monocultures (bootstrap mean [95%% CI], 1/day)\n');
    fprintf('  alpha untreated %.4f [%.4f, %.4f]   alpha +TGF %.4f [%.4f, %.4f]\n', ...
        M(1,1).rMean, M(1,1).rCI, M(1,2).rMean, M(1,2).rCI);
    fprintf('  delta untreated %.4f [%.4f, %.4f]   delta +TGF %.4f [%.4f, %.4f] (not used)\n', ...
        M(2,1).rMean, M(2,1).rCI, M(2,2).rMean, M(2,2).rCI);

    %% ---- co-cultures ----
    clear S
    for i = 1:size(raw,1)
        L = raw{i,3}(:);  N = raw{i,4}(:);
        L = L(1+s0:end);  N = N(1+s0:end);
        S(i).cond  = raw{i,1};                 %#ok<SAGROW>
        S(i).ratio = raw{i,2};
        S(i).x     = L ./ N;
        S(i).n     = N ./ N(1);
        S(i).alpha = alphaFix(S(i).cond);
        S(i).delta = deltaFix(S(i).cond);
    end
    tC   = (0:4-s0)';
    allX = cell2mat(arrayfun(@(s) s.x(2:end), S(:), 'UniformOutput', false));
    allN = cell2mat(arrayfun(@(s) s.n(2:end), S(:), 'UniformOutput', false));
    W    = [std(allX, 1), std(allN, 1)];
    nRes = 2 * (numel(tC) - 1) * numel(S);
    x0max = max(arrayfun(@(s) s.x(1), S));
    fprintf('Co-cultures: n = %d fitted points; largest initial leader fraction x0 = %.3f\n', nRes, x0max);

    aic = zeros(1,4);
    for h = 1:4
        [q, sse] = fit_case(h, S, tC, W, opts);
        kk = caseK(h) + 1;
        aic(h) = nRes*log(sse/nRes) + 2*kk + 2*kk*(kk+1)/(nRes - kk - 1);
        % residual bootstrap of this case
        c  = sqrt(nRes / (nRes - caseK(h)));
        PB = nan(nBoot, 4);
        for b = 1:nBoot
            try
                PB(b,:) = fit_case(h, boot_resample(S, q, tC, c), tC, W, opts);
            catch
            end
        end
        PB = PB(~any(isnan(PB), 2), :);
        R(h).q = q;  R(h).sse = sse;  R(h).ci = prctile(PB, [2.5 97.5], 1); %#ok<SAGROW>
        for cnd = 1:2
            a = alphaFix(cnd);  d = deltaFix(cnd);
            bq = q(2*cnd-1);    gq = q(2*cnd);
            R(h).xs(cnd)   = xstar(a, bq, gq, d);
            R(h).game{cnd} = game_name(a, bq, gq, d);
            xb = arrayfun(@(k) xstar(a, PB(k,2*cnd-1), PB(k,2*cnd), d), 1:size(PB,1));
            gb = arrayfun(@(k) game_name(a, PB(k,2*cnd-1), PB(k,2*cnd), d), 1:size(PB,1), 'UniformOutput', false);
            R(h).sameGame(cnd) = 100 * mean(strcmp(gb, R(h).game{cnd}));
            R(h).xsCI(cnd,:)   = prctile(xb, [2.5 97.5]);
        end
    end
    dAIC = aic - min(aic);
    [~, best] = min(aic);

    %% ---- print ----
    for h = 1:4
        fprintf('\n%s  (dAICc = %.2f)%s\n', caseNames{h}, dAIC(h), repmat('   <- best case', 1, h == best));
        fprintf('  beta  untreated %.3f [%.3f, %.3f]   +TGF %.3f [%.3f, %.3f]\n', ...
            R(h).q(1), R(h).ci(:,1), R(h).q(3), R(h).ci(:,3));
        fprintf('  gamma untreated %.3f [%.3f, %.3f]   +TGF %.3f [%.3f, %.3f]\n', ...
            R(h).q(2), R(h).ci(:,2), R(h).q(4), R(h).ci(:,4));
        for cnd = 1:2
            if any(strcmp(R(h).game{cnd}, {'Coordination', 'Snowdrift'}))
                xsTxt = sprintf('x* = %.2f [%.2f, %.2f]', R(h).xs(cnd), R(h).xsCI(cnd,:));
            else
                xsTxt = sprintf('x* = %.2f (outside [0,1]: no interior equilibrium)', R(h).xs(cnd));
            end
            fprintf('  %-10s %-19s %s; same class in %.1f%% of bootstrap fits\n', ...
                condTxt{cnd}, R(h).game{cnd}, xsTxt, R(h).sameGame(cnd));
        end
        fprintf('  stable outcome for the observed x0 (<= %.3f): %s\n', x0max, outcome_text(R(h), x0max));
    end

    OUT(cal).startDay = s0;  OUT(cal).alpha = alphaFix;  OUT(cal).delta = deltaFix; %#ok<SAGROW>
    OUT(cal).alphaCI = [M(1,1).rCI; M(1,2).rCI];  OUT(cal).deltaCI = M(2,1).rCI;
    OUT(cal).deltaTGF = [M(2,2).rMean M(2,2).rCI];
    OUT(cal).cases = R;  OUT(cal).dAIC = dAIC;  OUT(cal).best = best;  OUT(cal).x0max = x0max;
end
save('calibrated_params_both.mat', 'OUT');
fprintf('\nSaved calibrated_params_both.mat\n');

%% ================= Local functions =================
function R = mono_fit(T, t, nBoot, opts)
% Same as in the Fig5 script: y = y0*exp(r t), r in [0,1], y0 within the data range
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
    R.rMean = mean(P(:,1));
    R.rCI   = prctile(P(:,1), [2.5 97.5]);
end

function q = expand_params(p, h)
% [beta_U gamma_U beta_T gamma_T] for each case
    switch h
        case 1, q = [p(1) p(2) p(1) p(2)];      % H0: neither changes
        case 2, q = [p(1) p(2) p(3) p(2)];      % Hbeta: only beta changes
        case 3, q = [p(1) p(2) p(1) p(3)];      % Hgamma: only gamma changes
        otherwise, q = p(1:4);                  % Hbg: both change
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
        r = [r; (x(2:end) - S(i).x(2:end)) / W(1); (n(2:end) - S(i).n(2:end)) / W(2)]; %#ok<AGROW>
    end
end

function Sb = boot_resample(S, q, t, c)
% Residual bootstrap for single-replicate series (as in the Fig6 and FigS8 scripts)
    Sb = S;
    for i = 1:numel(S)
        [x, n] = lf_sim(q(2*S(i).cond-1 : 2*S(i).cond), S(i), t);
        rx = S(i).x(2:end) - x(2:end);  rx = (rx - mean(rx)) * c;
        rn = S(i).n(2:end) - n(2:end);  rn = (rn - mean(rn)) * c;
        m  = numel(rx);
        Sb(i).x = [S(i).x(1); x(2:end) + rx(randi(m, m, 1))];
        Sb(i).n = [1;         n(2:end) + rn(randi(m, m, 1))];
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
    [~, Z] = ode45(@(tt,z) rhs_2D(tt, z, p), tout, [Si.x(1); 1], odeset('RelTol',1e-8, 'AbsTol',1e-10));
    x = Z(:,1);  n = Z(:,2);
end

function xs = xstar(a, b, g, d)
    den = (b - d) + (g - a);
    if abs(den) < 1e-12, xs = NaN; else, xs = (b - d) / den; end
end

function g = game_name(a, b, gm, d)
    if b > d && gm > a,      g = 'Snowdrift';
    elseif b > d && gm < a,  g = 'Harmony';
    elseif b < d && gm < a,  g = 'Coordination';
    else,                    g = 'Prisoners Dilemma';
    end
end

function txt = outcome_text(Rh, x0max)
% Stable state reached from the observed initial leader fractions, per condition
    parts = cell(1, 2);
    for cnd = 1:2
        switch Rh.game{cnd}
            case 'Prisoners Dilemma', parts{cnd} = 'x -> 0 (followers only, from any x0)';
            case 'Harmony',           parts{cnd} = 'x -> 1 (leaders only, from any x0)';
            case 'Snowdrift',         parts{cnd} = sprintf('x -> x* = %.2f (coexistence)', Rh.xs(cnd));
            case 'Coordination'
                if x0max < Rh.xs(cnd)
                    parts{cnd} = sprintf('x -> 0 (followers only: all x0 < x* = %.2f)', Rh.xs(cnd));
                else
                    parts{cnd} = sprintf('depends on x0 (x* = %.2f)', Rh.xs(cnd));
                end
        end
    end
    txt = sprintf('untreated: %s | +TGF: %s', parts{1}, parts{2});
end
