%% calib_ferreira_emax.m
% Calibration of the Ferreira generator against the twin (twin_numerico):
% grid over KT (Tmax scaling) and Emax (ESPVR slope), LAP fixed.
% Question: does contractility (Emax) change the SHAPE of the ejection
% (shape factor k = SV/(Qpeak*LVET), time to peak), or only its scale?
% Targets from Charlton 2019 (PWDB baseline 25 y). Emin stays at the paper
% value (EDPVR; its effect is equivalent to LAP).
%
% REQUIRED in the InitFcn of twin_numerico, right after parameters_ferreira
% (use your variable names for Tmax, Emax, LAP):
%   if ~exist('KT', 'var'), KT = 1; end            % Tmax scaling (1 = Ferreira 2005)
%   Tmax = KT*Tmax;
%   if exist('LAP_SET', 'var'),  LAP  = LAP_SET;  end
%   if exist('EMAX_SET', 'var'), Emax = EMAX_SET; end
%
% Run from the Command Window: calib_ferreira_emax

MODEL  = 'twin_numerico';
T_STOP = '12';
LAP_FIX = 7.00;                  % [mmHg]
TGT = struct('SV', 67.1, 'Qpeak', 421, 'LVET', 280, 'tpk', 77, 'k', 0.57);   % Charlton
W   = [1 0.3 1];                 % cost weights: SV, Qpeak, LVET (as decided)

G4.KT   = 1.10:0.02:1.30;
G4.EMAX = 1.4:0.2:2.6;           % [mmHg/mL]
nK = numel(G4.KT);  nE = numel(G4.EMAX);  N = nK*nE;
G4.R = nan(nK, nE, 7);           % SV Qpeak tpeak LVET SBP DBP k
iSim = 0;  t0 = tic;               % (not "n": the InitFcn scripts clear n, k, s, r)
assignin('base', 'LAP_SET', LAP_FIX);
for a = 1:nK
    for b = 1:nE
        iSim = iSim + 1;
        assignin('base', 'KT', G4.KT(a));
        assignin('base', 'EMAX_SET', G4.EMAX(b));
        try
            out = sim(MODEL, 'StopTime', T_STOP);
            m = beat_metrics(out);
            G4.R(a,b,:) = [m, m(1)/(m(2)*1e-3*m(4))];
        catch err
            warning('KT %.2f, Emax %.2f failed: %s', G4.KT(a), G4.EMAX(b), err.message);
        end
        el = toc(t0);
        fprintf('>> [%d/%d, ~%.0f min left] KT %.2f  Emax %.2f : SV %5.1f  Qpeak %4.0f  tpk %3.0f  LVET %3.0f  P %5.1f/%5.1f  k %.2f\n', ...
                iSim, N, el/iSim*(N - iSim)/60, G4.KT(a), G4.EMAX(b), G4.R(a,b,1), G4.R(a,b,2), G4.R(a,b,3), ...
                G4.R(a,b,4), G4.R(a,b,5), G4.R(a,b,6), G4.R(a,b,7));
    end
end
evalin('base', 'clear KT LAP_SET EMAX_SET');
save('calib_ferreira_emax.mat', 'G4', 'TGT', 'LAP_FIX');

%% Best point for the chosen weights (re-run this section to change W)
E = W(1)*((G4.R(:,:,1) - TGT.SV)   /TGT.SV   ).^2 + ...
    W(2)*((G4.R(:,:,2) - TGT.Qpeak)/TGT.Qpeak).^2 + ...
    W(3)*((G4.R(:,:,4) - TGT.LVET) /TGT.LVET ).^2;
[~, m] = min(E(:));  [a, b] = ind2sub(size(E), m);
fprintf('\nBest (W = [%g %g %g], LAP %.2f): KT %.2f, Emax %.2f\n', W, LAP_FIX, G4.KT(a), G4.EMAX(b));
fprintf('  SV %.1f (%.1f)  Qpeak %.0f (%.0f)  tpk %.0f (%.0f)  LVET %.0f (%.0f)  P %.1f/%.1f  k %.2f (%.2f)\n', ...
        G4.R(a,b,1), TGT.SV, G4.R(a,b,2), TGT.Qpeak, G4.R(a,b,3), TGT.tpk, G4.R(a,b,4), TGT.LVET, ...
        G4.R(a,b,5), G4.R(a,b,6), G4.R(a,b,7), TGT.k);

%% Maps: black line = Charlton target
names = {'SV [mL]', 'Qpeak [mL/s]', 't peak [ms]', 'LVET [ms]', 'SBP [mmHg]', 'DBP [mmHg]', 'k (shape)'};
tg = [TGT.SV TGT.Qpeak TGT.tpk TGT.LVET NaN NaN TGT.k];
figure('Name', sprintf('Ferreira calibration - KT x Emax (LAP %.1f)', LAP_FIX));
for n = 1:7
    subplot(2, 4, n);
    imagesc(G4.EMAX, G4.KT, G4.R(:,:,n));  axis xy;  colorbar;
    xlabel('Emax [mmHg/mL]');  ylabel('KT');  title(names{n});
    if ~isnan(tg(n))
        hold on;  contour(G4.EMAX, G4.KT, G4.R(:,:,n), [tg(n) tg(n)], 'k', 'LineWidth', 2);
    end
end

%% ---------------------------------------------------------------------
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
