%% run_charlton_1d.m
% Offline run of the full PWDB network (Charlton 2019, baseline 25 y) with
% the explicit fixed-step 1D solver and the PWDB aortic inflow.
% Purpose: check that the solver reproduces the PWDB haemodynamics before
% using it in the twin.
%
% Needs in the folder: parameters_charlton.m, pwdb_geo0001.csv,
% make_charlton_inflow.m, AorticFlowWave_orig.m, charlton1d_setup.m,
% charlton1d_step.m
%
% Output in the workspace: res (time series at the measurement sites),
% S (network), A Q Pn Pc (final state).

%% ---------------- Settings ----------------
DT     = 1e-3;      % [s] time step. 1 ms: 15 short connectors (circle of Willis,
                    %     hand arches) are merged and 2 short terminals removed.
                    %     2e-4 s with DX_MAX = 0.01 keeps all 116 segments
                    %     (results differ by < 1 mmHg / < 0.3 mL/s).
DX_MAX = 0.02;      % [m] maximum cell length
T_END  = 10;        % [s] simulated time
T_REC  = 1e-3;      % [s] recording interval

%% ---------------- Parameters and inflow ----------------
parameters_charlton;
make_charlton_inflow;
prm = struct('rho', rho_ch, 'mu', mu_ch, 'alpha', alpha_ch, 'k1', k1_ch, 'k2', k2_ch, 'k3', k3_ch, ...
             'b0', b0_ch, 'b1', b1_ch, 'Pref', CH_DBP, 'Pout', P_venous);
S = charlton1d_setup(GEO_FILE, prm, DT, DX_MAX);

%% ---------------- Measurement sites (as in the PWDB) ----------------
% name, segment, position along the segment (0 = inlet, 1 = outlet)
sites = {'AorticRoot', 1, 0;  'Carotid', 15, 0.5;  'Brachial', 21, 0.5; ...
         'Radial', 22, 1;     'Femoral', 46, 0.5;  'AbdAorta', 39, 0.5; ...
         'SubclavianR_origin', 4, 0;  'CarotidR_origin', 5, 0; ...
         'CarotidL_origin', 15, 0;   'SubclavianL_origin', 19, 0};
nS = size(sites, 1);
gf = zeros(nS, 1); gc = gf;                 % face (flow) and cell (pressure) index
for k = 1:nS
    s = sites{k,2};
    gf(k) = S.fFirst(s) + round(sites{k,3}*(S.fLast(s) - S.fFirst(s)));
    gc(k) = S.cFirst(s) + round(sites{k,3}*(S.cLast(s) - S.cFirst(s)));
end

%% ---------------- Initial state: rest at P = DBP ----------------
A = S.A0;  Q = zeros(S.nF,1);  Pn = S.Pref*ones(S.nNodes,1);  Pc = S.Pref*ones(numel(S.tSeg),1);

%% ---------------- Time loop ----------------
nSteps = round(T_END/DT);  kRec = round(T_REC/DT);
nRec = floor(nSteps/kRec);
res.t = zeros(nRec,1); res.P = zeros(nRec,nS); res.Q = zeros(nRec,nS);
res.sites = sites(:,1)';
cflMax = 0; j = 0; tic;
for n = 1:nSteps
    tn  = n*DT;
    Qin = 1e-6*interp1(CH_INFLOW_t, CH_INFLOW_Q, mod(tn, CH_T), 'linear');
    [A, Q, Pn, Pc] = charlton1d_step(S, A, Q, Pn, Pc, Qin);
    if mod(n, kRec) == 0
        j = j + 1;
        Pg = S.Pref + S.beta(gc)./S.A0(gc).*(sqrt(A(gc)) - S.sqA0(gc));
        res.t(j) = tn;
        res.P(j,:) = Pg'/133.322;        % [mmHg] (cell nearest to the site)
        res.Q(j,:) = 1e6*Q(gf)';      % [mL/s] (face nearest to the site)
        u = (Q(S.cLF) + Q(S.cRF))/2./A;
        cflMax = max(cflMax, max((abs(u) + S.cK.*A.^0.25)*DT./S.dxc));
        if any(~isfinite(Q)), error('Solution diverged at t = %.3f s.', tn); end
    end
end
fprintf('Simulated %.1f s in %.0f s (%.1fx real time). Max CFL = %.2f\n', ...
        T_END, toc, toc/T_END, cflMax);

%% ---------------- Last beat: summary ----------------
nb = floor(res.t(end)/CH_T);                     % last complete beat
i  = res.t > (nb-1)*CH_T & res.t <= nb*CH_T;
tb = res.t(i);
fprintf('\n%-20s %7s %7s %7s %8s %8s %8s %9s\n', 'Site', 'SBP', 'DBP', 'MAP', ...
        'Qmean', 'Qmin', 'Qmax', 'back [%]');
for k = 1:nS
    p = res.P(i,k); q = res.Q(i,k);
    fprintf('%-20s %7.1f %7.1f %7.1f %8.2f %8.2f %8.2f %9.0f\n', sites{k,1}, ...
            max(p), min(p), trapz(tb,p)/(tb(end)-tb(1)), trapz(tb,q)/(tb(end)-tb(1)), ...
            min(q), max(q), 100*max(0,-min(q))/max(q));
end
fprintf('PWDB baseline 25 y: DBP %.1f, MAP %.1f mmHg, CO %.2f L/min\n', ...
        CH_DBP, CH_MBP, CH_SV*CH_HR/1000);

%% ---------------- Plots: last beat at the PWDB sites ----------------
k6 = 1:6;
figure('Name', '1D PWDB network: last beat');
subplot(1,2,1); plot(tb - tb(1), res.P(i,k6)); grid on
xlabel('t [s]'); ylabel('P [mmHg]'); legend(sites(k6,1), 'Interpreter', 'none');
subplot(1,2,2); plot(tb - tb(1), res.Q(i,k6)); grid on
xlabel('t [s]'); ylabel('Q [mL/s]'); yline(0, 'k:');
