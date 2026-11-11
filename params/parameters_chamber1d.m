%% parameters_chamber1d.m
% Builds the 1D sub-tree of each Windkessel chamber (PWDB network, Charlton
% 2019) for the MATLAB Function blocks chamber1d_upper / _abdominal / _lower.
% Run AFTER parameters_charlton.m (e.g. in init_model).
%
% Each block: input Q measured into the chamber [mL/s], output P_target [mmHg].
% Block sample time = DT_1D (fixed step of the 1D solver).

DT_1D     = 1e-3;     % [s] time step of the 1D solver (= sample time of the blocks)
DX_MAX_1D = 0.02;     % [m] maximum cell length

prm_1d = struct('rho', rho_ch, 'mu', mu_ch, 'alpha', alpha_ch, 'k1', k1_ch, 'k2', k2_ch, ...
                'k3', k3_ch, 'b0', b0_ch, 'b1', b1_ch, 'Pref', CH_DBP, 'Pout', P_venous);
stubF_1d = [F_STUB_SUBCLAVIAN_R F_STUB_CAROTID_R F_STUB_CAROTID_L F_STUB_SUBCLAVIAN_L];

S_upper     = charlton1d_setup(GEO_FILE, prm_1d, DT_1D, DX_MAX_1D, 'upper', stubF_1d);
S_abdominal = charlton1d_setup(GEO_FILE, prm_1d, DT_1D, DX_MAX_1D, 'abdominal');
S_lower     = charlton1d_setup(GEO_FILE, prm_1d, DT_1D, DX_MAX_1D, 'lower');

%% Full network (twin: aorta of Charlton + all trees), for the block "TwinFull1D"
% Input: LVP through the aortic valve. Outputs: aortic-root pressure, and the
% flow into / pressure at each chamber, at the same places as the bench
% (after the phantom stubs). Used as reference for the bench.
% Aortic stiffness of the reference twin (the aorta is the phantom on the
% bench): K_AORTA multiplies Eh of the aortic segments (1 = Charlton).
% K_AORTA = 2 -> aorta half as compliant. Does NOT change the chamber blocks.
AORTA_SEGS = [1 2 3 14 18 27 28 35 37 39 41];   % asc., arch A/B, innominate, thoracic A/B, abdominal A-E
K_AORTA    = 1;
prm_full = prm_1d;  prm_full.Kstiff = ones(size(CH_R(:)));  prm_full.Kstiff(AORTA_SEGS) = K_AORTA;
S_full = charlton1d_setup(GEO_FILE, prm_full, DT_1D, DX_MAX_1D, 'full');
rootSegs = {[4 5 15 19], [26 29 34 36 38 40], [42 43]};
rootF    = {stubF_1d, zeros(1,6), zeros(1,2)};
nm = {'Up', 'Ab', 'Lo'};
for r = 1:3                             % faces (flow) and cells (pressure) at the chamber inlets
    s = rootSegs{r};  k = round(rootF{r}(:).*S_full.N(s(:)));
    S_full.(['f' nm{r}]) = S_full.fFirst(s(:)) + k;
    S_full.(['c' nm{r}]) = S_full.cFirst(s(:)) + min(k, S_full.N(s(:)) - 1);
end
clear rootSegs rootF nm r s k

R_V = 0.003;                       % [mmHg*s/mL] aortic valve resistance (one-way valve on the bench)
S_full.Rv = R_V*133.322e6;         % -> [Pa*s/m^3]
