%% calib_ferreira_grid.m
% Calibration of the Ferreira generator against the twin (twin_numerico):
% grid over KT (Tmax scaling) and LAP, targets from Charlton 2019 (PWDB
% baseline 25 y): SV, peak aortic flow, LVET.
%
% Two stages:
%   1) coarse grid over the whole range (KT 1.00:0.05:1.35, LAP 7:1:12)
%   2) fine grid around the best coarse point (KT step 0.01, LAP step 0.25)
% Every point is simulated once and saved, so the weights W can be changed
% afterwards (section "Best point") without re-running.
%
% REQUIRED in the InitFcn of twin_numerico, right after parameters_ferreira:
%   if ~exist('KT', 'var'), KT = 1; end           % Tmax scaling (1 = Ferreira 2005)
%   Tmax = KT*Tmax;
%   if exist('LAP_SET', 'var'), LAP = LAP_SET; end % LAP override (calibration only)
%
% Run from the Command Window: calib_ferreira_grid

MODEL  = 'twin_numerico';
T_STOP = '12';                   % [s] simulated time per point (transient settled)
TGT = struct('SV', 67.1, 'Qpeak', 421, 'LVET', 280);   % Charlton: mL, mL/s, ms
W   = [1 0.3 0.6];                   % cost weights: SV, Qpeak, LVET
names = {'SV [mL]', 'Qpeak [mL/s]', 't peak [ms]', 'LVET [ms]', 'SBP [mmHg]', 'DBP [mmHg]'};

%% Stage 1: coarse grid
G1.KT  = 1.00:0.05:1.35;
G1.LAP = 7:1:12;
G1.R   = run_grid(MODEL, T_STOP, G1.KT, G1.LAP);
[a, b] = best_point(G1.R, TGT, W);
save('calib_ferreira_grid.mat', 'G1', 'TGT', 'names');

%% Stage 2: fine grid around the best coarse point
G2.KT  = max(1.00, G1.KT(a) - 0.05) : 0.01 : G1.KT(a) + 0.05;
G2.LAP = max(5, G1.LAP(b) - 1) : 0.25 : min(12, G1.LAP(b) + 1);
G2.R   = run_grid(MODEL, T_STOP, G2.KT, G2.LAP);
save('calib_ferreira_grid.mat', 'G1', 'G2', 'TGT', 'names');

%% Best point for the chosen weights (re-run this section to change W)
[a, b] = best_point(G2.R, TGT, W);
fprintf('\nBest (W = [%g %g %g]): KT %.2f, LAP %.2f\n', W, G2.KT(a), G2.LAP(b));
fprintf('  SV %.1f (%.1f)  Qpeak %.0f (%.0f)  LVET %.0f (%.0f)  P %.1f/%.1f\n', ...
        G2.R(a,b,1), TGT.SV, G2.R(a,b,2), TGT.Qpeak, G2.R(a,b,4), TGT.LVET, G2.R(a,b,5), G2.R(a,b,6));

%% Maps (coarse and fine): black line = Charlton target
tg = [TGT.SV TGT.Qpeak NaN TGT.LVET NaN NaN];
G = {G1, G2};  ttl = {'coarse', 'fine'};
for g = 1:2
    figure('Name', ['Ferreira calibration - ' ttl{g}]);
    for n = 1:6
        subplot(2, 3, n);
        imagesc(G{g}.LAP, G{g}.KT, G{g}.R(:,:,n));  axis xy;  colorbar;
        xlabel('LAP [mmHg]');  ylabel('KT');  title(names{n});
        if ~isnan(tg(n))
            hold on;  contour(G{g}.LAP, G{g}.KT, G{g}.R(:,:,n), [tg(n) tg(n)], 'k', 'LineWidth', 2);
        end
    end
end
evalin('base', 'clear KT LAP_SET');

%% ---------------------------------------------------------------------
function R = run_grid(model, tStop, ktGrid, lapGrid)
% Simulates every (KT, LAP) pair; KT and LAP_SET go to the base workspace,
% where the InitFcn runs.
nK = numel(ktGrid);  nL = numel(lapGrid);  N = nK*nL;
R  = nan(nK, nL, 6);  n = 0;  t0 = tic;
for a = 1:nK
    for b = 1:nL
        n = n + 1;
        assignin('base', 'KT', ktGrid(a));
        assignin('base', 'LAP_SET', lapGrid(b));
        try
            out = sim(model, 'StopTime', tStop);
            R(a,b,:) = beat_metrics(out);
        catch err
            warning('KT %.2f, LAP %.2f failed: %s', ktGrid(a), lapGrid(b), err.message);
        end
        el = toc(t0);
        fprintf('>> [%d/%d, ~%.0f min left] KT %.2f  LAP %5.2f : SV %5.1f  Qpeak %4.0f  tpk %3.0f  LVET %3.0f  P %5.1f/%5.1f\n', ...
                n, N, el/n*(N - n)/60, ktGrid(a), lapGrid(b), R(a,b,1), R(a,b,2), R(a,b,3), R(a,b,4), R(a,b,5), R(a,b,6));
    end
end
end

function [a, b] = best_point(R, TGT, W)
E = W(1)*((R(:,:,1) - TGT.SV)   /TGT.SV   ).^2 + ...
    W(2)*((R(:,:,2) - TGT.Qpeak)/TGT.Qpeak).^2 + ...
    W(3)*((R(:,:,4) - TGT.LVET) /TGT.LVET ).^2;
[~, m] = min(E(:));  [a, b] = ind2sub(size(E), m);
end

function m = beat_metrics(out)
% Last complete beat, from one valve opening to the next.
Qs = out.logsout.get('Q_ao').Values;  Ps = out.logsout.get('P_root').Values;
q  = squeeze(Qs.Data);  t = Qs.Time;  p = squeeze(Ps.Data);
op = find(diff(q > 1) == 1) + 1;
i  = op(end-1):op(end)-1;
qt = q(i);  tt = t(i);  pt = p(i);
[qmax, k] = max(qt);  c = find(qt > 1, 1, 'last');
m = [trapz(tt, qt), qmax, 1e3*(tt(k) - tt(1)), 1e3*(tt(c) - tt(1)), max(pt), min(pt)];
end
