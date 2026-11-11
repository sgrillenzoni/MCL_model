%% test_chamber_1d.m
% Offline check of the three chamber sub-trees (upper, abdominal, lower)
% BEFORE putting them in Simulink.
%
% 1) Runs the full 1D network with the PWDB inflow and records, for each
%    chamber, the flow that enters its sub-tree and the pressure there.
% 2) Drives each sub-tree ALONE with that recorded flow (as the measured
%    flow will drive the chamber block on the bench).
% 3) Compares the sub-tree root pressure (= chamber pressure target) with
%    the pressure of the full network at the same place.
% Lower: exact check (one inlet node). Abdominal and upper: the inlets sit
% at different aortic points and are joined into one node (the chamber), so
% small differences are expected.
%
% Needs: parameters_charlton.m, pwdb_geo0001.csv, make_charlton_inflow.m,
% AorticFlowWave_orig.m, charlton1d_setup.m, charlton1d_step.m

DT = 1e-3;  DX_MAX = 0.02;  T_END = 8;

parameters_charlton;
make_charlton_inflow;
prm = struct('rho', rho_ch, 'mu', mu_ch, 'alpha', alpha_ch, 'k1', k1_ch, 'k2', k2_ch, ...
             'k3', k3_ch, 'b0', b0_ch, 'b1', b1_ch, 'Pref', CH_DBP, 'Pout', P_venous);
stubF = [F_STUB_SUBCLAVIAN_R F_STUB_CAROTID_R F_STUB_CAROTID_L F_STUB_SUBCLAVIAN_L];

%% ---------------- 1) full network, record chamber inflows ----------------
Sf = charlton1d_setup(GEO_FILE, prm, DT, DX_MAX, 'full');
reg = {'upper', [4 5 15 19], stubF;  'abdominal', [26 29 34 36 38 40], zeros(1,6); ...
       'lower', [42 43], zeros(1,2)};
for r = 1:3                                   % face and cell at the inlet of each root
    segs = reg{r,2};  Fr = reg{r,3};
    reg{r,4} = Sf.fFirst(segs) + round(Fr(:).*Sf.N(segs));                    % faces
    reg{r,5} = Sf.cFirst(segs) + min(round(Fr(:).*Sf.N(segs)), Sf.N(segs)-1); % cells
end
nSteps = round(T_END/DT);
Qrec = zeros(nSteps, 3);  Prec = zeros(nSteps, 3);
A = Sf.A0; Q = zeros(Sf.nF,1); Pn = Sf.Pref*ones(Sf.nNodes,1); Pc = Sf.Pref*ones(numel(Sf.tSeg),1);
tic
for n = 1:nSteps
    Qin = 1e-6*interp1(CH_INFLOW_t, CH_INFLOW_Q, mod(n*DT, CH_T), 'linear');
    [A, Q, Pn, Pc] = charlton1d_step(Sf, A, Q, Pn, Pc, Qin);
    for r = 1:3
        Qrec(n,r) = sum(Q(reg{r,4}));
        ic = reg{r,5};
        Prec(n,r) = mean(Sf.Pref + Sf.beta(ic)./Sf.A0(ic).*(sqrt(A(ic)) - Sf.sqA0(ic)));
    end
end
fprintf('Full network: %.1f s simulated in %.1f s\n', T_END, toc);

%% ---------------- 2) each sub-tree alone, driven by the recorded flow ----------------
t  = (1:nSteps)'*DT;
nb = floor(T_END/CH_T);  i = t > (nb-1)*CH_T & t <= nb*CH_T;      % last beat
fprintf('\n%-10s %-28s %-28s %s\n', 'Chamber', 'P sub-tree SBP/DBP/MAP', 'P full net SBP/DBP/MAP', 'RMSE [mmHg]');
Pch = zeros(nSteps, 3);
for r = 1:3
    Sr = charlton1d_setup(GEO_FILE, prm, DT, DX_MAX, reg{r,1}, reg{r,3});
    A = Sr.A0; Q = zeros(Sr.nF,1); Pn = Sr.Pref*ones(Sr.nNodes,1); Pc = Sr.Pref*ones(numel(Sr.tSeg),1);
    tic
    for n = 1:nSteps
        [A, Q, Pn, Pc] = charlton1d_step(Sr, A, Q, Pn, Pc, Qrec(n,r));
        Pch(n,r) = Pn(Sr.root);
    end
    tr = toc;
    p1 = Pch(i,r)/133.322;  p2 = Prec(i,r)/133.322;
    fprintf('%-10s %6.1f / %5.1f / %5.1f         %6.1f / %5.1f / %5.1f         %5.2f   (%.1f s)\n', ...
            reg{r,1}, max(p1), min(p1), mean(p1), max(p2), min(p2), mean(p2), ...
            sqrt(mean((p1 - p2).^2)), tr);
end

figure('Name', 'Chamber sub-trees vs full network');
for r = 1:3
    subplot(1,3,r);
    plot(t(i) - t(find(i,1)), Prec(i,r)/133.322, 'k', t(i) - t(find(i,1)), Pch(i,r)/133.322, 'r--');
    grid on; title(reg{r,1}); xlabel('t [s]');
    if r == 1, ylabel('P [mmHg]'); legend('full network', 'sub-tree (chamber target)'); end
end

%% ---------------- signals for the Simulink test (step 5) ----------------
% Recorded inflow [mL/s] and full-network pressure [mmHg] of each chamber,
% as timeseries for "From Workspace" blocks.
REF_upper     = timeseries([1e6*Qrec(:,1) Prec(:,1)/133.322], t);
REF_abdominal = timeseries([1e6*Qrec(:,2) Prec(:,2)/133.322], t);
REF_lower     = timeseries([1e6*Qrec(:,3) Prec(:,3)/133.322], t);
fprintf('\nREF_upper, REF_abdominal, REF_lower created: column 1 = Q [mL/s], column 2 = P [mmHg].\n');
