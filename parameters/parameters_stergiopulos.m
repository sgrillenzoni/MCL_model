%% parameters_stergiopulos.m
% Lumped-parameter values for the arterial tree of Stergiopulos et al. (1992).
%
% Source:
%   Stergiopulos N, Young DF, Rogge TR. Computer simulation of arterial flow
%   with applications to arterial and aortic stenoses. J Biomech 1992;25(12):1477-1488.
%   Table 1 (segment geometry and compliance) and Table 2 (terminal impedances).
%
% Electrical analogy: voltage = pressure [mmHg], current = flow [mL/s].
% Units:
%   R, RT  [mmHg*s/mL]      resistance
%   L      [mmHg*s^2/mL]    inertance
%   C, CT  [mL/mmHg]        compliance
%
% Naming:  <type>_<vessel>[_<side>]
%   R_, L_, C_   arterial segments (55)
%   RT_, CT_     terminal (peripheral) loads of the terminal segments (28),
%                modelled as a 2-element Windkessel at the distal end
%   Letters A/B/C... mark consecutive sub-segments of the same vessel.
%
% How the values were obtained (CGS units, then converted to mmHg and mL):
%   - Resistance (Poiseuille):     R = 8*mu*l / (pi*r^4)
%   - Inertance:                   L = rho*l / (pi*r^2)
%     For tapered segments (proximal radius r1 different from distal r2) the
%     radius is assumed to vary linearly with the axial position x along the
%     segment, from r1 at the inlet (x = 0) to r2 at the outlet (x = l):
%                                  r(x) = r1 + (r2 - r1) * x / l
%     The segment is treated as a series of infinitesimal cylindrical slices
%     of length dx, each obeying the formulas above. Summing the slices in
%     series gives R = 8*mu/pi * int_0^l dx/r(x)^4 and
%     L = rho/pi * int_0^l dx/r(x)^2, which evaluate to:
%                                  R = 8*mu/pi * l*(r1^2 + r1*r2 + r2^2) / (3*r1^3*r2^3)
%                                  L = rho/pi  * l / (r1*r2)
%     (both reduce to the expressions above when r1 = r2).
%     Conversion: R[mmHg*s/mL] = R[dyn*s/cm^5] / 1333.22,
%                 L[mmHg*s^2/mL] = L[g/cm^4] / 1333.22.
%   - Blood properties: viscosity mu = 0.04 P (0.004 Pa*s), density
%     rho = 1.05 g/cm^3. These are intentionally blood values: the model
%     represents the arterial tree with blood, independently of the working
%     fluid used in the physical circuit.
%   - C taken directly from Table 1 (volumetric compliance), unit-converted only.
%   - RT and CT converted from the SI values of Table 2
%     (1 mmHg = 133.322 Pa, 1 mL = 1e-6 m^3).
%   - Segment 20 is labelled only 'Vertebral' in the source table; it branches
%     from the left subclavian artery, so it is named vertebral_left here.
%
% Usage: run this script before the simulation, or call it from the model's
% PreLoadFcn / InitFcn callback (Model Properties > Callbacks). In each
% Simscape block, enter the variable name (e.g. R_renal_left) as the value.

%% ===================== Arterial segments ======================

%% Segments 1-17: ascending aorta, aortic arch, right arm, carotids

% Ascending aorta (segment 1)
R_ascendingAorta = 6.821e-05;
L_ascendingAorta = 0.0004737;
C_ascendingAorta = 0.1392;

% Aortic arch A (segment 2)
R_aorticArchA = 9.711e-05;
L_aorticArchA = 0.0003997;
C_aorticArchA = 0.03946;

% Innominate (segment 3)
R_innominate = 0.001758;
L_innominate = 0.002217;
C_innominate = 0.018;

% R. subclavian A (segment 4)
R_subclavianA_right = 0.008114;
L_subclavianA_right = 0.004764;
C_subclavianA_right = 0.007466;

% R. carotid (segment 5)
R_commonCarotid_right = 0.07215;
L_commonCarotid_right = 0.03241;
C_commonCarotid_right = 0.02848;

% R. vertebral (segment 6)  - terminal
R_vertebral_right = 0.9555;
L_vertebral_right = 0.1078;
C_vertebral_right = 0.002242;

% R. subclavian B (segment 7)
R_subclavianB_right = 0.3913;
L_subclavianB_right = 0.1112;
C_subclavianB_right = 0.04516;

% R. radial (segment 8)  - terminal
R_radial_right = 2.982;
L_radial_right = 0.2384;
C_radial_right = 0.002502;

% R. ulnar A (segment 9)
R_ulnarA_right = 0.2396;
L_ulnarA_right = 0.03634;
C_ulnarA_right = 0.00148;

% R. interosseous (segment 10)  - terminal
R_interosseous_right = 8.802;
L_interosseous_right = 0.2392;
C_interosseous_right = 0.00012;

% R. ulnar B (segment 11)  - terminal
R_ulnarB_right = 0.9501;
L_ulnarB_right = 0.1154;
C_ulnarB_right = 0.002946;

% R. internal carotid (segment 12)  - terminal
R_internalCarotid_right = 7.522;
L_internalCarotid_right = 0.302;
C_internalCarotid_right = 0.001257;

% R. external carotid (segment 13)  - terminal
R_externalCarotid_right = 7.522;
L_externalCarotid_right = 0.302;
C_externalCarotid_right = 0.001257;

% Aortic arch B (segment 14)
R_aorticArchB = 0.0002273;
L_aorticArchB = 0.000854;
C_aorticArchB = 0.06946;

% L. carotid (segment 15)
R_commonCarotid_left = 0.08479;
L_commonCarotid_left = 0.03809;
C_commonCarotid_left = 0.03346;

% L. internal carotid (segment 16)  - terminal
R_internalCarotid_left = 7.522;
L_internalCarotid_left = 0.302;
C_internalCarotid_left = 0.001257;

% L. external carotid (segment 17)  - terminal
R_externalCarotid_left = 7.522;
L_externalCarotid_left = 0.302;
C_externalCarotid_left = 0.001257;

%% Segments 18-27: thoracic aorta, left arm, intercostals

% Thoracic aorta A (segment 18)
R_thoracicAortaA = 0.0003989;
L_thoracicAortaA = 0.001306;
C_thoracicAortaA = 0.07959;

% L. subclavian A (segment 19)
R_subclavianA_left = 0.008114;
L_subclavianA_left = 0.004764;
C_subclavianA_left = 0.007466;

% Vertebral (segment 20)  - terminal
R_vertebral_left = 0.9555;
L_vertebral_left = 0.1078;
C_vertebral_left = 0.002242;

% L. subclavian B (segment 21)
R_subclavianB_left = 0.3913;
L_subclavianB_left = 0.1112;
C_subclavianB_left = 0.04516;

% L. radial (segment 22)  - terminal
R_radial_left = 2.982;
L_radial_left = 0.2384;
C_radial_left = 0.002502;

% L. ulnar A (segment 23)
R_ulnarA_left = 0.2396;
L_ulnarA_left = 0.03634;
C_ulnarA_left = 0.00148;

% L. interosseous (segment 24)  - terminal
R_interosseous_left = 8.802;
L_interosseous_left = 0.2392;
C_interosseous_left = 0.00012;

% L. ulnar B (segment 25)  - terminal
R_ulnarB_left = 0.9501;
L_ulnarB_left = 0.1154;
C_ulnarB_left = 0.002946;

% Intercostals (segment 26)  - terminal
R_intercostals = 0.698;
L_intercostals = 0.06685;
C_intercostals = 0.004;

% Thoracic aorta B (segment 27)
R_thoracicAortaB = 0.004195;
L_thoracicAortaB = 0.005988;
C_thoracicAortaB = 0.06346;

%% Segments 28-41: abdominal aorta, visceral and renal branches

% Abdominal aorta A (segment 28)
R_abdominalAortaA = 0.002925;
L_abdominalAortaA = 0.003571;
C_abdominalAortaA = 0.0272;

% Celiac A (segment 29)
R_celiacA = 0.003302;
L_celiacA = 0.001648;
C_celiacA = 0.001813;

% Celiac B (segment 30)
R_celiacB = 0.04775;
L_celiacB = 0.006267;
C_celiacB = 0.001333;

% Hepatic (segment 31)  - terminal
R_hepatic = 0.2153;
L_hepatic = 0.03419;
C_hepatic = 0.003066;

% Gastric (segment 32)  - terminal
R_gastric = 0.5167;
L_gastric = 0.05494;
C_gastric = 0.002013;

% Splenic (segment 33)  - terminal
R_splenic = 0.08416;
L_splenic = 0.02088;
C_splenic = 0.004986;

% Superior mesenteric (segment 34)  - terminal
R_superiorMesenteric = 0.01259;
L_superiorMesenteric = 0.007816;
C_superiorMesenteric = 0.01387;

% Abdominal aorta B (segment 35)
R_abdominalAortaB = 0.0005895;
L_abdominalAortaB = 0.0006964;
C_abdominalAortaB = 0.005333;

% L. renal (segment 36)  - terminal
R_renal_left = 0.0535;
L_renal_left = 0.01187;
C_renal_left = 0.002226;

% Abdominal aorta C (segment 37)
R_abdominalAortaC = 0.0006305;
L_abdominalAortaC = 0.0007202;
C_abdominalAortaC = 0.005066;

% R. renal (segment 38)  - terminal
R_renal_right = 0.0535;
L_renal_right = 0.01187;
C_renal_right = 0.002226;

% Abdominal aorta D (segment 39)
R_abdominalAortaD = 0.008025;
L_abdominalAortaD = 0.008361;
C_abdominalAortaD = 0.0452;

% Inferior mesenteric (segment 40)  - terminal
R_inferiorMesenteric = 0.5829;
L_inferiorMesenteric = 0.04896;
C_inferiorMesenteric = 0.001056;

% Abdominal aorta E (segment 41)
R_abdominalAortaE = 0.001045;
L_abdominalAortaE = 0.0009271;
C_abdominalAortaE = 0.004666;

%% Segments 42-55: iliac arteries and legs

% R. common iliac (segment 42)
R_commonIliac_right = 0.02673;
L_commonIliac_right = 0.01129;
C_commonIliac_right = 0.006106;

% L. common iliac (segment 43)
R_commonIliac_left = 0.02673;
L_commonIliac_left = 0.01129;
C_commonIliac_left = 0.006106;

% L. external iliac (segment 44)
R_externalIliac_left = 0.1488;
L_externalIliac_left = 0.04178;
C_externalIliac_left = 0.02082;

% L. internal iliac (segment 45)  - terminal
R_internalIliac_left = 0.2388;
L_internalIliac_left = 0.03134;
C_internalIliac_left = 0.0044;

% L. femoral (segment 46)
R_femoral_left = 1.443;
L_femoral_left = 0.2257;
C_femoral_left = 0.01819;

% L. deep femoral (segment 47)  - terminal
R_deepFemoral_left = 0.4422;
L_deepFemoral_left = 0.0666;
C_deepFemoral_left = 0.001507;

% L. posterior tibial (segment 48)  - terminal
R_posteriorTibial_left = 2.239;
L_posteriorTibial_left = 0.2311;
C_posteriorTibial_left = 0.002941;

% L. anterior tibial (segment 49)  - terminal
R_anteriorTibial_left = 9.175;
L_anteriorTibial_left = 0.5088;
C_anteriorTibial_left = 0.001123;

% R. external iliac (segment 50)
R_externalIliac_right = 0.1488;
L_externalIliac_right = 0.04178;
C_externalIliac_right = 0.02082;

% R. internal iliac (segment 51)  - terminal
R_internalIliac_right = 0.2388;
L_internalIliac_right = 0.03134;
C_internalIliac_right = 0.0044;

% R. femoral (segment 52)
R_femoral_right = 1.443;
L_femoral_right = 0.2257;
C_femoral_right = 0.01819;

% R. deep femoral (segment 53)  - terminal
R_deepFemoral_right = 0.4422;
L_deepFemoral_right = 0.0666;
C_deepFemoral_right = 0.001507;

% R. posterior tibial (segment 54)  - terminal
R_posteriorTibial_right = 2.239;
L_posteriorTibial_right = 0.2311;
C_posteriorTibial_right = 0.002941;

% R. anterior tibial (segment 55)  - terminal
R_anteriorTibial_right = 9.175;
L_anteriorTibial_right = 0.5088;
C_anteriorTibial_right = 0.001123;

%% ================= Terminal (peripheral) loads ==================
% Resistance RT from the distal node to the venous pressure node,
% compliance CT from the distal node to ground.
% Note: Stergiopulos chose CT proportional to 1/RT, so every terminal bed
% has the same time constant RT*CT = 0.186 s.

% R. vertebral (segment 6)
RT_vertebral_right = 45.08;
CT_vertebral_right = 0.004127;

% R. radial (segment 8)
RT_radial_right = 39.6;
CT_radial_right = 0.004698;

% R. interosseous (segment 10)
RT_interosseous_right = 632.3;
CT_interosseous_right = 0.0002942;

% R. ulnar B (segment 11)
RT_ulnarB_right = 39.6;
CT_ulnarB_right = 0.004698;

% R. internal carotid (segment 12)
RT_internalCarotid_right = 104.3;
CT_internalCarotid_right = 0.001784;

% R. external carotid (segment 13)
RT_externalCarotid_right = 104.3;
CT_externalCarotid_right = 0.001784;

% L. internal carotid (segment 16)
RT_internalCarotid_left = 104.3;
CT_internalCarotid_left = 0.001784;

% L. external carotid (segment 17)
RT_externalCarotid_left = 104.3;
CT_externalCarotid_left = 0.001784;

% Vertebral (segment 20)
RT_vertebral_left = 45.08;
CT_vertebral_left = 0.004127;

% L. radial (segment 22)
RT_radial_left = 39.6;
CT_radial_left = 0.004698;

% L. interosseous (segment 24)
RT_interosseous_left = 632.3;
CT_interosseous_left = 0.0002942;

% L. ulnar B (segment 25)
RT_ulnarB_left = 39.6;
CT_ulnarB_left = 0.004698;

% Intercostals (segment 26)
RT_intercostals = 10.43;
CT_intercostals = 0.01784;

% Hepatic (segment 31)
RT_hepatic = 27.23;
CT_hepatic = 0.006833;

% Gastric (segment 32)
RT_gastric = 40.58;
CT_gastric = 0.004585;

% Splenic (segment 33)
RT_splenic = 17.4;
CT_splenic = 0.01069;

% Superior mesenteric (segment 34)
RT_superiorMesenteric = 6.976;
CT_superiorMesenteric = 0.02667;

% L. renal (segment 36)
RT_renal_left = 8.476;
CT_renal_left = 0.02195;

% R. renal (segment 38)
RT_renal_right = 8.476;
CT_renal_right = 0.02195;

% Inferior mesenteric (segment 40)
RT_inferiorMesenteric = 51.6;
CT_inferiorMesenteric = 0.003605;

% L. internal iliac (segment 45)
RT_internalIliac_left = 59.53;
CT_internalIliac_left = 0.003125;

% L. deep femoral (segment 47)
RT_deepFemoral_left = 35.78;
CT_deepFemoral_left = 0.0052;

% L. posterior tibial (segment 48)
RT_posteriorTibial_left = 35.78;
CT_posteriorTibial_left = 0.0052;

% L. anterior tibial (segment 49)
RT_anteriorTibial_left = 41.93;
CT_anteriorTibial_left = 0.004437;

% R. internal iliac (segment 51)
RT_internalIliac_right = 59.53;
CT_internalIliac_right = 0.003125;

% R. deep femoral (segment 53)
RT_deepFemoral_right = 35.78;
CT_deepFemoral_right = 0.0052;

% R. posterior tibial (segment 54)
RT_posteriorTibial_right = 35.78;
CT_posteriorTibial_right = 0.0052;

% R. anterior tibial (segment 55)
RT_anteriorTibial_right = 41.93;
CT_anteriorTibial_right = 0.004437;

%% ===================== Venous pressure ========================
P_venous = 8;   % [mmHg] Venous pressure at the outlet of all terminal resistances RT (downstream end of the peripheral beds; physiological range ~2-8 mmHg). Use a single source shared by all 28 terminals.

%% Total arterial compliance (segments + terminals)
% Lumped equivalent for a 2-element Windkessel: at low frequency all the
% compliances charge together, so the equivalent is their sum.
clear C_total                          % avoid summing a previous result
vC = [who('C_*'); who('CT_*')];
C_total = 0;
for k = 1:numel(vC)
    C_total = C_total + eval(vC{k});
end
clear vC k
fprintf('Total arterial compliance: %.3f mL/mmHg\n', C_total);
%% Sanity check
% All terminal resistances in parallel should give the total peripheral
% resistance, about 1 mmHg*s/mL.
RT_all = [ ...
    RT_vertebral_right, ...
    RT_radial_right, ...
    RT_interosseous_right, ...
    RT_ulnarB_right, ...
    RT_internalCarotid_right, ...
    RT_externalCarotid_right, ...
    RT_internalCarotid_left, ...
    RT_externalCarotid_left, ...
    RT_vertebral_left, ...
    RT_radial_left, ...
    RT_interosseous_left, ...
    RT_ulnarB_left, ...
    RT_intercostals, ...
    RT_hepatic, ...
    RT_gastric, ...
    RT_splenic, ...
    RT_superiorMesenteric, ...
    RT_renal_left, ...
    RT_renal_right, ...
    RT_inferiorMesenteric, ...
    RT_internalIliac_left, ...
    RT_deepFemoral_left, ...
    RT_posteriorTibial_left, ...
    RT_anteriorTibial_left, ...
    RT_internalIliac_right, ...
    RT_deepFemoral_right, ...
    RT_posteriorTibial_right, ...
    RT_anteriorTibial_right, ...
];
fprintf('Stergiopulos parameters loaded: 55 segments, %d terminal loads.\n', numel(RT_all));
fprintf('Total peripheral resistance (terminals in parallel): %.3f mmHg*s/mL\n', 1/sum(1./RT_all));
R_periph = 1/sum(1./RT_all);   % [mmHg*s/mL] total peripheral resistance
clear RT_all
