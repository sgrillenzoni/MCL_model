%% parameters_charlton.m
% Lumped-parameter values of the arterial network of the Pulse Wave Database
% (PWDB), baseline 25-year-old subject, converted for the Simscape twin.
%
% Source:
%   Charlton PH, Mariscal Harana J, Vennin S, Li Y, Chowienczyk P, Alastruey J.
%   Modeling arterial pulse waves in healthy aging: a database for in silico
%   evaluation of hemodynamics and pulse wave indexes.
%   Am J Physiol Heart Circ Physiol 2019;317(5):H1062-H1085.
%   Geometry: pwdb_geo0001.csv (116 segments). Global parameters:
%   pwdb_model_configs.csv, subject 1 (base = 1, age 25).
%   Formulas as in the PWDB code (create_input_files_for_a_simulation.m).
%
% Electrical analogy: voltage = pressure [mmHg], current = flow [mL/s].
% Units: R [mmHg*s/mL], L [mmHg*s^2/mL], C [mL/mmHg], TAU [s].
%
% How the values are obtained (SI, then converted):
%   radius varies linearly along each segment; the formulas below are
%   integrated along the segment (201 points, trapezoidal rule).
%   - Resistance (parabolic profile, alpha = 4/3):  R' = 8*mu/(pi*r^4)
%   - Inertance:                                     L' = rho/(pi*r^2)
%   - Wall stiffness (PWDB):  Eh/r = 0.1*(k1*exp(k2*100*r) + k3)   [Pa]
%     Compliance per length (thin wall, linearised):  C' = 1.5*pi*r^3/Eh
%     (equivalent to C' = A/(rho*c^2), with c^2 = 2/3*(Eh/r)/rho)
%   - Wall viscosity (PWDB, Mynard et al. 2015):
%       Gamma = b0 + b1/D   [g/s],  D = diameter [cm]
%     In each cell this is a Voigt element: Rv = TAU/C_cell, with
%       TAU = 3*Gamma/(2*sqrt(pi)*Eh)   (same for every cell of a segment)
%   - Terminal beds: 3-element Windkessel, as in the PWDB:
%       R1 = Z0 = rho*c/A at the segment outlet (characteristic impedance)
%       R2 = R_total - R1,  C = peripheral compliance
%     Outflow pressure p_out = 33.2 mmHg (PWDB value, used as P_venous).
%
% Note: the given radii are at diastolic pressure; the linearisation uses
% them directly (small difference from the reference radius of the 1D model).
%
% Output in the workspace:
%   CH_* vectors indexed by segment number (1..116), used by
%   build_charlton_tree.m; named R_/L_/C_/TAU_ variables for the aortic and
%   phantom segments, same names as parameters_stergiopulos.m, so the
%   hand-built aortic blocks keep working.
%
% Usage: run INSTEAD of parameters_stergiopulos.m, parameters_vascular_tuning.m
% and parameters_phantom.m (model callbacks: InitFcn / PreLoadFcn).

%% ---------------- Files and subject ----------------
GEO_FILE = 'pwdb_geo0001.csv';     % baseline 25-year-old subject

% pwdb_model_configs.csv, subject 1 (base = 1, age 25)
CH_HR   = 72.916;     % [bpm]
CH_SV   = 66.334;     % [mL]
CH_PFT  = 0.079;      % [s]   time of peak aortic flow
CH_LVET = 0.282;      % [s]   ejection time
CH_RFV  = 0.73;       % [mL]  reverse flow volume at the aortic root
CH_DBP  = 74.9;       % [mmHg]
CH_MBP  = 88.925;     % [mmHg]
mu_ch   = 0.0025;     % [Pa*s]  blood viscosity
alpha_ch = 4/3;      % [-]     velocity profile coefficient (PWDB: 1.3333)
rho_ch  = 1060;       % [kg/m^3]
k1_ch = 3e6;  k2_ch = -13.5;  k3_ch = 5.3849e5;   % stiffness [g/s^2/cm], [1/cm], [g/s^2/cm]
b0_ch = 600;  b1_ch = 150;                         % viscosity [g/s], [g*cm/s]

P_venous = 33.2;      % [mmHg] outflow pressure of all terminal Windkessels (PWDB p_out)

%% ---------------- Tuning factors (1 = PWDB values) ----------------
KTAU_V     = 1;       % scale of the wall viscosity time constants
KC_VIRTUAL = 1;       % scale of the compliances of the virtual trees
TAU_V = 0.004;        % [s] only for the hand-built aortic Rv block(s), see note at the end

%% ---------------- Read the geometry ----------------
G = csvread(GEO_FILE, 1, 0);       % skip header (csvread: available in R2018)
nSeg   = size(G, 1);
CH_IN  = G(:, 2)';   CH_OUT = G(:, 3)';
len_m  = G(:, 4)';   r_in = G(:, 5)';   r_out = G(:, 6)';
Cp_SI  = G(:, 7)';   Rp_SI = G(:, 8)';
CH_LEN  = 100*len_m;               % [cm], used for the number of cells
CH_TERM = Rp_SI > 0;               % terminal segments

%% ---------------- Segment properties ----------------
R_SI2MM = 133.322*1e6;             % 1 mmHg*s/mL = 1.33322e8 Pa*s/m^3
x  = linspace(0, 1, 201);
CH_R = zeros(1, nSeg); CH_L = CH_R; CH_C = CH_R; CH_TAU = CH_R;
for s = 1:nSeg
    r   = r_in(s) + (r_out(s) - r_in(s))*x;          % [m]
    A   = pi*r.^2;
    Ehr = 0.1*(k1_ch*exp(k2_ch*100*r) + k3_ch);      % Eh/r [Pa]
    Eh  = Ehr.*r;                                     % [Pa*m]
    Gam = (b0_ch + b1_ch./(200*r))/1000;             % [kg/s]
    CH_R(s)   = len_m(s)*trapz(x, 8*mu_ch./(pi*r.^4)) / R_SI2MM;
    CH_L(s)   = len_m(s)*trapz(x, rho_ch./A)          / R_SI2MM;
    CH_C(s)   = len_m(s)*trapz(x, 1.5*pi*r.^3./Eh)    * R_SI2MM;
    CH_TAU(s) = mean(3*Gam./(2*sqrt(pi)*Eh));         % [s]
end

%% ---------------- Terminal Windkessels ----------------
c_out  = sqrt(2/3/rho_ch*0.1*(k1_ch*exp(k2_ch*100*r_out) + k3_ch));   % [m/s]
Z0     = rho_ch*c_out./(pi*r_out.^2) / R_SI2MM;                         % [mmHg*s/mL]
CH_R1  = zeros(1, nSeg); CH_R2 = CH_R1; CH_CT = CH_R1;
CH_R1(CH_TERM) = Z0(CH_TERM);
CH_R2(CH_TERM) = Rp_SI(CH_TERM)/R_SI2MM - Z0(CH_TERM);
CH_CT(CH_TERM) = Cp_SI(CH_TERM)*R_SI2MM;
if any(CH_R2(CH_TERM) <= 0)
    error('parameters_charlton: R2 <= 0 for segments %s.', mat2str(find(CH_TERM & CH_R2 <= 0)));
end

%% ---------------- Named variables (aorta and phantom segments) ----------------
% Same names as parameters_stergiopulos.m, for the hand-built aortic blocks
% and parameters_phantom.m. Numbering as in the PWDB (= Stergiopulos 1-55).
named = {1,'ascendingAorta'; 2,'aorticArchA'; 3,'innominate'; 4,'subclavianA_right'; ...
         5,'commonCarotid_right'; 14,'aorticArchB'; 15,'commonCarotid_left'; ...
         18,'thoracicAortaA'; 19,'subclavianA_left'; 27,'thoracicAortaB'; ...
         28,'abdominalAortaA'; 35,'abdominalAortaB'; 37,'abdominalAortaC'; ...
         39,'abdominalAortaD'; 41,'abdominalAortaE'};
for k = 1:size(named, 1)
    s = named{k,1};  n = named{k,2};
    eval(sprintf('R_%s = CH_R(%d); L_%s = CH_L(%d); C_%s = CH_C(%d); TAU_%s = CH_TAU(%d);', ...
                 n, s, n, s, n, s, n, s));
end

%% ---------------- Phantom stubs (physical) and virtual remainders ----------------
% The phantom has four supra-aortic outlets: R subclavian and R common
% carotid (after the innominate, which is fully physical), L common carotid,
% L subclavian. The first centimetres of each ("stub") are physical; the rest
% of the segment is virtual and is the ROOT of the upper-body tree.
% Replaces parameters_phantom.m (same variable names, PWDB lengths).
% !!! STUB LENGTHS ARE PLACEHOLDERS: measure them on the phantom
% !!! (from the arch wall to the barbed connector).
L_STUB_SUBCLAVIAN_R = 2.0;   % [cm] seg 4,  R subclavian A   (PWDB length 3.4 cm)
L_STUB_CAROTID_R    = 3.0;   % [cm] seg 5,  R common carotid (PWDB length 9.4 cm)
L_STUB_CAROTID_L    = 3.0;   % [cm] seg 15, L common carotid (PWDB length 13.9 cm)
L_STUB_SUBCLAVIAN_L = 2.0;   % [cm] seg 19, L subclavian A   (PWDB length 3.4 cm),
                             %      outlet BEFORE the L vertebral origin
% Physical fraction of each segment, limited to [0, 1]
F_STUB_SUBCLAVIAN_R = min(max(L_STUB_SUBCLAVIAN_R/CH_LEN(4),  0), 1);
F_STUB_CAROTID_R    = min(max(L_STUB_CAROTID_R   /CH_LEN(5),  0), 1);
F_STUB_CAROTID_L    = min(max(L_STUB_CAROTID_L   /CH_LEN(15), 0), 1);
F_STUB_SUBCLAVIAN_L = min(max(L_STUB_SUBCLAVIAN_L/CH_LEN(19), 0), 1);
% Stub (outside the block) and remainder (inside the block), as in
% parameters_phantom.m: R_stub_<v>, R_rest_<v>, L_..., C_...
stub = {'subclavian_R', 4, F_STUB_SUBCLAVIAN_R; 'carotid_R', 5, F_STUB_CAROTID_R; ...
        'carotid_L', 15, F_STUB_CAROTID_L; 'subclavian_L', 19, F_STUB_SUBCLAVIAN_L};
for k = 1:size(stub, 1)
    v = stub{k,1};  s = stub{k,2};  f = stub{k,3};
    eval(sprintf('R_stub_%s = %.10g*CH_R(%d); R_rest_%s = %.10g*CH_R(%d);', v, f, s, v, 1-f, s));
    eval(sprintf('L_stub_%s = %.10g*CH_L(%d); L_rest_%s = %.10g*CH_L(%d);', v, f, s, v, 1-f, s));
    eval(sprintf('C_stub_%s = %.10g*CH_C(%d); C_rest_%s = %.10g*CH_C(%d);', v, f, s, v, 1-f, s));
end
fprintf('Phantom stubs (physical fraction): subclavian R %.2f, carotid R %.2f, carotid L %.2f, subclavian L %.2f\n', ...
        F_STUB_SUBCLAVIAN_R, F_STUB_CAROTID_R, F_STUB_CAROTID_L, F_STUB_SUBCLAVIAN_L);

%% ---------------- Globals for parameters_ferreira.m ----------------
% Used ONLY by the original Ferreira block (internal arterial load R1, C3,
% for its standalone verification); the twin's heart_chamber does not use them.
R_periph = 1/sum(1./(CH_R1(CH_TERM) + CH_R2(CH_TERM)));   % [mmHg*s/mL] terminals in parallel
C_total  = sum(CH_C) + sum(CH_CT);                         % [mL/mmHg] vessels + terminals

%% ---------------- Summary ----------------
fprintf('PWDB (Charlton 2019) baseline 25 y loaded: %d segments, %d terminal Windkessels.\n', ...
        nSeg, nnz(CH_TERM));
fprintf('Total peripheral resistance: %.3f mmHg*s/mL   (PWDB pvr: 0.691)\n', ...
        1/sum(1./(CH_R1(CH_TERM) + CH_R2(CH_TERM))));
fprintf('Compliance: vessels %.2f + terminals %.2f = %.2f mL/mmHg\n', ...
        sum(CH_C), sum(CH_CT), sum(CH_C) + sum(CH_CT));
fprintf('Expected mean flow (MBP - p_out)/TPR = %.1f mL/s   (SV*HR/60 = %.1f mL/s)\n', ...
        (CH_MBP - P_venous)*sum(1./(CH_R1(CH_TERM) + CH_R2(CH_TERM))), CH_SV*CH_HR/60);

% Note on TAU_V: the hand-built ascending-aorta cell uses
% Rv = max(TAU_V,1e-9)/C_ascendingAorta. The PWDB value for that segment is
% TAU_ascendingAorta (~0.6 ms); TAU_V = 4 ms is kept there only to damp the
% 20 Hz ringing seen with the Ferreira ventricle. Write TAU_ascendingAorta in
% that block to use the PWDB value.

clear G x r A Ehr Eh Gam s k n v f stub named len_m r_in r_out Cp_SI Rp_SI c_out Z0 R_SI2MM nSeg
