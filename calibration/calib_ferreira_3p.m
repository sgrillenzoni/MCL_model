%% calib_ferreira_3p.m
% Final calibration of the Ferreira generator against the twin (twin_numerico).
% Each parameter has its own role:
%   Emax (ESPVR, contractility) -> shape of the ejection (k = SV/(Qpeak*LVET))
%   KT   (Tmax scaling)         -> duration of the ejection (LVET)
%   LAP  (filling)              -> stroke volume (scales SV and peak together)
% Targets from Charlton 2019 (PWDB baseline 25 y); weights fixed beforehand.
%
% REQUIRED in the InitFcn of twin_numerico, right after parameters_ferreira:
%   if ~exist('KT', 'var'), KT = 1; end
%   Tmax = KT*Tmax;
%   if exist('LAP_SET', 'var'),  LAP  = LAP_SET;  end
%   if exist('EMAX_SET', 'var'), Emax = EMAX_SET; end
%
% NOTE: this script runs in the base workspace, like the InitFcn, whose
% scripts clear several short names (n, k, s, r, v, f, x, A, G...).
% Loop variables here are iE, iK, iL, iSim: do not rename them.
%
% Run from the Command Window: calib_ferreira_3p

MODEL  = 'twin_numerico';
T_STOP = '12';
TGT = struct('SV', 67.1, 'Qpeak', 421, 'LVET', 280, 'tpk', 77, 'k', 0.57);   % Charlton
W   = [1 0.3 1];                 % cost weights: SV, Qpeak, LVET (as decided)

G5.EMAX = [1.4 1.5 1.6 1.8];         % [mmHg/mL]
G5.KT   = 1.12:0.01:1.22;
G5.LAP  = 7:0.5:10;              % [mmHg]
nE = numel(G5.EMAX);  nK = numel(G5.KT);  nL = numel(G5.LAP);  N = nE*nK*nL;
G5.R = nan(nE, nK, nL, 7);       % SV Qpeak tpeak LVET SBP DBP k
iSim = 0;  t0 = tic;
for iE = 1:nE
    for iK = 1:nK
        for iL = 1:nL
            iSim = iSim + 1;
            assignin('base', 'EMAX_SET', G5.EMAX(iE));
            assignin('base', 'KT',       G5.KT(iK));
            assignin('base', 'LAP_SET',  G5.LAP(iL));
            try
                out = sim(MODEL, 'StopTime', T_STOP);
                mm  = beat_metrics(out);
                G5.R(iE,iK,iL,:) = [mm, mm(1)/(mm(2)*1e-3*mm(4))];
            catch err
                warning('Emax %.2f, KT %.2f, LAP %.2f failed: %s', G5.EMAX(iE), G5.KT(iK), G5.LAP(iL), err.message);
            end
            el = toc(t0);
            fprintf('>> [%d/%d, ~%.0f min left] Emax %.2f KT %.2f LAP %5.2f : SV %5.1f Qpeak %4.0f tpk %3.0f LVET %3.0f P %5.1f/%5.1f k %.2f\n', ...
                    iSim, N, el/iSim*(N - iSim)/60, G5.EMAX(iE), G5.KT(iK), G5.LAP(iL), squeeze(G5.R(iE,iK,iL,:)));
        end
    end
end
evalin('base', 'clear KT LAP_SET EMAX_SET');
save('calib_ferreira_3p.mat', 'G5', 'TGT');

%% Best point for the chosen weights (re-run this section to change W)
C = W(1)*((G5.R(:,:,:,1) - TGT.SV)   /TGT.SV   ).^2 + ...
    W(2)*((G5.R(:,:,:,2) - TGT.Qpeak)/TGT.Qpeak).^2 + ...
    W(3)*((G5.R(:,:,:,4) - TGT.LVET) /TGT.LVET ).^2;
fprintf('\nW = [%g %g %g]\n', W);
for iE = 1:nE                    % best for each Emax, then overall
    Ce = squeeze(C(iE,:,:));  [~, mi] = min(Ce(:));  [iK, iL] = ind2sub(size(Ce), mi);
    rr = squeeze(G5.R(iE,iK,iL,:));
    fprintf('Emax %.2f: KT %.2f, LAP %.2f -> SV %.1f  Qpeak %.0f  tpk %.0f  LVET %.0f  P %.1f/%.1f  k %.2f  (cost %.4f)\n', ...
            G5.EMAX(iE), G5.KT(iK), G5.LAP(iL), rr, Ce(mi));
end
[~, mi] = min(C(:));  [iE, iK, iL] = ind2sub(size(C), mi);
fprintf('BEST: Emax %.2f, KT %.2f, LAP %.2f\n', G5.EMAX(iE), G5.KT(iK), G5.LAP(iL));
fprintf('Charlton: SV %.1f  Qpeak %.0f  tpk %.0f  LVET %.0f  P 101.7/76.6  k %.2f\n', ...
        TGT.SV, TGT.Qpeak, TGT.tpk, TGT.LVET, TGT.k);

%% Maps for each Emax (KT x LAP): black line = Charlton target
names = {'SV [mL]', 'Qpeak [mL/s]', 't peak [ms]', 'LVET [ms]', 'SBP [mmHg]', 'DBP [mmHg]', 'k (shape)'};
tg = [TGT.SV TGT.Qpeak TGT.tpk TGT.LVET NaN NaN TGT.k];
for iE = 1:nE
    figure('Name', sprintf('Ferreira calibration - Emax %.2f', G5.EMAX(iE)));
    for iM = 1:7
        subplot(2, 4, iM);
        Z = squeeze(G5.R(iE,:,:,iM));
        imagesc(G5.LAP, G5.KT, Z);  axis xy;  colorbar;
        xlabel('LAP [mmHg]');  ylabel('KT');  title(sprintf('%s, Emax %.1f', names{iM}, G5.EMAX(iE)));
        if ~isnan(tg(iM))
            hold on;  contour(G5.LAP, G5.KT, Z, [tg(iM) tg(iM)], 'k', 'LineWidth', 2);
        end
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
