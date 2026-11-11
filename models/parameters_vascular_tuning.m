%% parameters_vascular_tuning.m
% Tuning factors of the virtual arterial trees built by build_vascular_tree.m.
% Run AFTER parameters_stergiopulos.m. All factors = original network when set
% to the "original" values below; change them here, no need to rebuild blocks.
%
% Effects found with the offline replica of the twin (Python, Nov 2026),
% compared by eye with Charlton et al. 2019 (25-year-old waveforms):
%   - cells <= 3 cm (build option) removes the artificial resonances of long
%     vessels; mean flows and flow split unchanged
%   - TAU_V 3-5 ms: carotid backflow 55% -> ~39% of peak, abdominal flow
%     always forward, less ringing; mean flows unchanged
%   - KCT_LEGS = 3: femoral backflow ~40% -> ~17% of peak (Charlton ~15%)
%   - KRT_HEAD < 1: less carotid backflow, but more flow to the upper body
%     (0.25 gives 33% of CO to the upper body: too much; ~0.5 should be a
%     compromise, NOT tested)
% TAU_V, KCT_LEGS and KRT_HEAD are MODIFICATIONS of the Stergiopulos network:
% justify the chosen values in the report.

%% Wall viscosity (Voigt model): Rv = TAU_V / C for every cell [s]
TAU_V = 0.004;          % original: 0 (purely elastic). Tested: 0.003-0.005
                    % !!! value to be justified with literature on arterial wall viscoelasticity

%% Compliance scale of the VIRTUAL vessels and terminals [-]
KC_VIRTUAL = 1;     % original: 1. (Aortic = physical compliance: from the phantom measurement)

%% Cerebral terminal beds (vertebrals, internal and external carotids) [-]
% RT scaled by KRT_HEAD, CT by 1/KRT_HEAD (time constant RT*CT unchanged)
KRT_HEAD = 1;       % original: 1. Tested: 0.6, 0.35, 0.25
%% Upper body terminal beds: resistance factors [-]
% RT scaled by k, CT by 1/k (time constant RT*CT unchanged)
KRT_CAROTID = 0.2;    % internal + external carotids. original: 1
KRT_VERT    = 1;      % vertebrals. original: 1
KRT_ARMS    = 2;      % radial, ulnar B, interosseous. original: 1

%% Arm terminal beds: compliance factor (RT unchanged) [-]
KCT_ARMS = 6;   % original: 1
for v = {'radial_right','ulnarB_right','interosseous_right', ...
         'radial_left','ulnarB_left','interosseous_left'}
    eval(sprintf('CT_%s = KCT_ARMS*CT_%s;', v{1}, v{1}));
end

grp = {KRT_CAROTID, {'internalCarotid_right','externalCarotid_right','internalCarotid_left','externalCarotid_left'}; ...
       KRT_VERT,    {'vertebral_right','vertebral_left'}; ...
       KRT_ARMS,    {'radial_right','ulnarB_right','interosseous_right', ...
                     'radial_left','ulnarB_left','interosseous_left'}};
for g = 1:size(grp, 1)
    for v = grp{g, 2}
        eval(sprintf('RT_%s = %g*RT_%s;  CT_%s = CT_%s/%g;', ...
                     v{1}, grp{g,1}, v{1}, v{1}, v{1}, grp{g,1}));
    end
end
%% Leg terminal beds: CT scaled by KCT_LEGS (RT unchanged) [-]
KCT_LEGS = 1;       % original: 1. Tested: 3

%% Left common carotid: fraction that is physical (stub in the phantom) [-]
F_STUB_CAROTID_L = 0.15;   % !!! placeholder: measure the stub on the phantom

fprintf('Vascular tuning: TAU_V = %g s, KC_VIRTUAL = %g, KRT_HEAD = %g, KCT_LEGS = %g, F_STUB_CAROTID_L = %g\n', ...
        TAU_V, KC_VIRTUAL, KRT_HEAD, KCT_LEGS, F_STUB_CAROTID_L);
%% TEST: scale of the aortic compliances [-]
KC_AORTA = 1;   % original: 1
aorta = {'ascendingAorta','aorticArchA','aorticArchB','thoracicAortaA','thoracicAortaB', ...
         'abdominalAortaA','abdominalAortaB','abdominalAortaC','abdominalAortaD','abdominalAortaE'};
for k = 1:numel(aorta)
    eval(sprintf('C_%s = KC_AORTA*C_%s;', aorta{k}, aorta{k}));
end