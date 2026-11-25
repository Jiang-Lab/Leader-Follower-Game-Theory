%% ========================================================================
%  Supplementary Fig. S10: SENSITIVITY + REGIME ROBUSTNESS  (label fig:sens1)
%    A-D : local OAT (+/-10%) sensitivity of x(t) and N(t), both conditions
%    E-F : regime robustness in the (beta-delta, gamma-alpha) plane
%  Growth-coupled 2D model.  Run once per K (set K_cells below).
%
%  Parameters p = [alpha beta gamma delta] come from the calibration in
%  the Fig5 script (calibrated_params.mat):
%    alpha  : leader monoculture rate, one per condition (TGF-beta1 lowers it)
%    delta  : follower monoculture rate, one value shared by both conditions
%    beta, gamma : shared by both conditions (case H0, joint fit of all ratios)
%  The calibration data stay below 1e7 cells, so N/K < 1e-2 for K = 1e9 and
%  < 1e-3 for K = 1e10: the same parameters apply for both K values.
%% ========================================================================
clear; clc; close all;
rng(1,'twister'); warning('off','all'); format short;

% -------------------- USER SETTING: choose ONE ---------------------------
%K_cells = 1e9;              % <-- 1e9 for one run, 1e10 for the other
K_cells = 1e10;
ratioIC = 1;                 % co-culture used for A-D: 1 = 1:1, 2 = 1:2, 3 = 1:5 L:F
% -------------------------------------------------------------------------

pLocal  = 0.10;    % +/-10% OAT
pGlobal = 0.25;    % +/-25% joint Monte Carlo
Nmc     = 4000;
tdata   = [0;1;2;3;4];

% -------- calibrated parameters (from the Fig5 script) ----------------
if isfile('calibrated_params.mat')
    cal = load('calibrated_params.mat');
    alphaFix = cal.alphaFix;  deltaFix = cal.deltaFix;  bg = cal.bgHat;
    fprintf('Loaded calibrated_params.mat\n');
else
    % values from the joint fit reported in the paper (update if refitted)
    alphaFix = [0.3792 0.1212];      % untreated, +TGF-beta1
    deltaFix = [0.4470 0.4470];      % shared
    bg       = [0.123  0.766];       % [beta gamma], shared
    fprintf('calibrated_params.mat not found: using the values typed in the script\n');
end
C(1) = struct('name','Without TGF-\beta1','p',[alphaFix(1) bg(1) bg(2) deltaFix(1)]);
C(2) = struct('name','With TGF-\beta1',   'p',[alphaFix(2) bg(1) bg(2) deltaFix(2)]);

% co-culture initial conditions (x0 = L0/N0, N0); same seeding in both conditions
L0 = [4.5e5 3.6e5 1.8e5];   N0 = [8.7e5 1.0e6 9.9e5];       % 1:1, 1:2, 1:5 L:F
ratioTxt = {'1:1','1:2','1:5'};
for c = 1:2
    C(c).x0 = L0(ratioIC)/N0(ratioIC);  C(c).N0 = N0(ratioIC);
end
fprintf('Initial condition for A-D: %s L:F co-culture (x0 = %.3f)\n', ratioTxt{ratioIC}, C(1).x0);
for c = 1:2
    fprintf('%-18s alpha = %.4f  beta = %.4f  gamma = %.4f  delta = %.4f\n', ...
            strrep(C(c).name,'\',''), C(c).p);
end

% -------------------- style (as in the other figures) --------------------
baseFont = 'Helvetica';
LABEL_FS = 18;  TITLE_FS = 16;  TICK_FS = 11;  LEGEND_FS = 12;  PANEL_FS = 20;  NOTE_FS = 12;
set(groot,'DefaultAxesFontName',baseFont,'DefaultTextFontName',baseFont, ...
    'DefaultAxesTickDir','out','DefaultAxesBox','off','DefaultAxesLineWidth',1);

param_names = {'\alpha','\beta','\gamma','\delta'};
colDown=[0.00 0.45 0.74]; colUp=[0.85 0.33 0.10];
colPD=[.565 .933 .565]; colCD=[.678 .847 .902]; colH=[.855 .439 .839]; colSD=[.957 .647 .510];
opts = odeset('NonNegative',1:2,'RelTol',1e-8,'AbsTol',1e-10);

fig = figure('Color','w','Units','inches','Position',[1 1 12.5 15]);
panelLab = {'A','B','C','D','E','F'};

% ======================= PART A-D : local OAT ============================
for c = 1:2
    p0 = C(c).p;  z0 = [C(c).x0; C(c).N0];
    [~,Zb] = ode23s(@(t,z) rhs2d(t,z,p0,K_cells), tdata, z0, opts);
    xb = Zb(:,1); Nb = Zb(:,2);

    SIx = zeros(4,2); SIN = zeros(4,2);
    for i = 1:4
        for s = 1:2
            pp = p0; pp(i) = p0(i)*(1 + (2*s-3)*pLocal);   % -,+
            [~,Z] = ode23s(@(t,z) rhs2d(t,z,pp,K_cells), tdata, z0, opts);
            SIx(i,s) = norm(Z(:,1)-xb)/max(norm(xb),eps);
            SIN(i,s) = norm(Z(:,2)-Nb)/max(norm(Nb),eps);
        end
    end

    % leader-fraction panel
    ax = subplot(3,2,2*c-1); b=bar(SIx,'grouped');
    b(1).FaceColor=colDown; b(2).FaceColor=colUp; b(1).EdgeColor='none'; b(2).EdgeColor='none';
    grid on; ax.GridAlpha = 0.18;
    set(ax,'XTickLabel',param_names,'FontSize',TICK_FS);
    ylabel('Sensitivity Index','FontWeight','bold','FontSize',LABEL_FS); ylim([0 max(SIx(:))*1.35]);
    title(sprintf('Leader Fraction x(t), %s',C(c).name),'FontWeight','bold','FontSize',TITLE_FS);
    legend({'-10%','+10%'},'Location','northwest','Box','off','FontSize',LEGEND_FS);
    panelLetter(ax,panelLab{2*c-1},PANEL_FS);

    % population panel
    ax = subplot(3,2,2*c); b=bar(SIN,'grouped');
    b(1).FaceColor=colDown; b(2).FaceColor=colUp; b(1).EdgeColor='none'; b(2).EdgeColor='none';
    grid on; ax.GridAlpha = 0.18;
    set(ax,'XTickLabel',param_names,'FontSize',TICK_FS);
    ylabel('Sensitivity Index','FontWeight','bold','FontSize',LABEL_FS); ylim([0 max(SIN(:))*1.35]);
    title(sprintf('Cell Population N(t), %s',C(c).name),'FontWeight','bold','FontSize',TITLE_FS);
    panelLetter(ax,panelLab{2*c},PANEL_FS);

    [~,ix]=sort(max(SIx,[],2),'descend'); [~,iN]=sort(max(SIN,[],2),'descend');
    fprintf('%s: ranking  x(t): %s  |  N(t): %s\n', strrep(C(c).name,'\',''), ...
        strjoin(strrep(param_names(ix),'\',''),' > '), strjoin(strrep(param_names(iN),'\',''),' > '));
end

% ======================= PART E-F : regime plane =========================
for c = 1:2
    p0=C(c).p; a0=p0(1); b0=p0(2); g0=p0(3); d0=p0(4);
    ax=subplot(3,2,4+c); hold(ax,'on'); L=1.05;
    % regions (colors pre-blended with white instead of transparency, for fast export)
    rectangle('Position',[-L 0 L L],'FaceColor',colPD*0.55+0.45,'EdgeColor','none');
    rectangle('Position',[0 0 L L],'FaceColor',colSD*0.55+0.45,'EdgeColor','none');
    rectangle('Position',[-L -L L L],'FaceColor',colCD*0.55+0.45,'EdgeColor','none');
    rectangle('Position',[0 -L L L],'FaceColor',colH*0.55+0.45,'EdgeColor','none');
    text(-0.72,0.92,sprintf('Prisoner''s\nDilemma'),'HorizontalAlignment','center','FontSize',NOTE_FS,'FontWeight','bold','Color',[.25 .25 .25]);
    text(0.62,0.92,'Snowdrift','HorizontalAlignment','center','FontSize',NOTE_FS,'FontWeight','bold','Color',[.25 .25 .25]);
    text(-0.72,-0.85,'Coordination','HorizontalAlignment','center','FontSize',NOTE_FS,'FontWeight','bold','Color',[.25 .25 .25]);
    text(0.62,-0.85,'Harmony','HorizontalAlignment','center','FontSize',NOTE_FS,'FontWeight','bold','Color',[.25 .25 .25]);

    mults=(1-pGlobal)+2*pGlobal*rand(Nmc,4); P=p0.*mults;
    plot(P(:,2)-P(:,4),P(:,3)-P(:,1),'.','Color',[.55 .55 .55],'MarkerSize',3);   % opaque: fast vector export
    fprintf('%s: +/-%d%% Monte Carlo, fraction still in Prisoner''s Dilemma = %.3f\n', ...
        strrep(C(c).name,'\',''), round(100*pGlobal), mean(P(:,2) < P(:,4) & P(:,3) > P(:,1)));

    for i=1:4, for s=[1-pLocal 1+pLocal]
        q=p0; q(i)=p0(i)*s;
        plot(q(2)-q(4),q(3)-q(1),'o','MarkerSize',6,'MarkerFaceColor',colUp,'MarkerEdgeColor','k','LineWidth',0.5);
    end, end
    plot(b0-d0,g0-a0,'p','MarkerSize',20,'MarkerFaceColor',[1 .84 0],'MarkerEdgeColor','k','LineWidth',1.2);
    plot([b0-d0 0],[g0-a0 g0-a0],'k-','LineWidth',1.3);
    plot([b0-d0 b0-d0],[g0-a0 0],'k-','LineWidth',1.3);
    text((b0-d0)/2,g0-a0+0.07,sprintf('%.2f',abs(b0-d0)),'HorizontalAlignment','center','FontSize',NOTE_FS,'FontWeight','bold');
    text(b0-d0-0.06,(g0-a0)/2,sprintf('%.2f',abs(g0-a0)),'HorizontalAlignment','right','FontSize',NOTE_FS,'FontWeight','bold','Rotation',90);
    plot([-L L],[0 0],'k-','LineWidth',1.6); plot([0 0],[-L L],'k-','LineWidth',1.6);
    xlim([-L L]); ylim([-L L]);
    xlabel('\beta - \delta','FontWeight','bold','FontSize',LABEL_FS);
    ylabel('\gamma - \alpha','FontWeight','bold','FontSize',LABEL_FS);
    set(ax,'FontSize',TICK_FS,'Layer','top'); title(C(c).name,'FontWeight','bold','FontSize',TITLE_FS);
    panelLetter(ax,panelLab{4+c},PANEL_FS);
    % legend below each regime panel (outside the axes, so it never covers the data)
    hFit = findobj(ax,'Marker','pentagram');
    hOAT = findobj(ax,'Marker','o','-and','MarkerFaceColor',colUp);
    hMC  = findobj(ax,'Marker','.');
    legend([hFit(1); hOAT(1); hMC(1)], ...
           {'Fitted', '\pm10% OAT', sprintf('\\pm%d%% Monte Carlo', round(100*pGlobal))}, ...
           'Location','southoutside','Orientation','horizontal','Box','off','FontSize',LEGEND_FS);

    sym_bd=(d0-b0)/(d0+b0)*100; sym_ga=(g0-a0)/(g0+a0)*100;
    if sym_bd<sym_ga, bind='beta/delta'; else, bind='gamma/alpha'; end
    fprintf('%s | K=%.0e : boundary  beta=delta %.1f%%  gamma=alpha %.1f%%  BINDING %s %.1f%%\n', ...
        strrep(C(c).name,'\',''),K_cells,sym_bd,sym_ga,bind,min(sym_bd,sym_ga));
end

drawnow;                                    % show the finished figure before saving
fprintf('Saving figure...\n');
outDir = fullfile('..', 'figures');  if ~isfolder(outDir), mkdir(outDir); end
exportgraphics(fig,fullfile(outDir,sprintf('FigS10_sensitivity_K1e%d.png',round(log10(K_cells)))),'Resolution',300);
exportgraphics(fig,fullfile(outDir,sprintf('FigS10_sensitivity_K1e%d.pdf',round(log10(K_cells)))),'ContentType','vector');
fprintf('\nSaved to %s\n',outDir);

%% ========================= LOCAL FUNCTIONS ==============================
function dzdt = rhs2d(~, z, p, K)
    x=z(1); N=z(2);
    A=p(1)*x+p(2)*(1-x); B=p(3)*x+p(4)*(1-x); h=1-N/K;   % N/K logistic
    dzdt = [ (A-B)*x*(1-x)*h ; (A*x+B*(1-x))*N*h ];
end
function panelLetter(ax,L,fs)
    text(ax,-0.16,1.12,L,'Units','normalized','FontSize',fs,'FontWeight','bold', ...
        'HorizontalAlignment','left','VerticalAlignment','top');
end
