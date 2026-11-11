%% ferreira_params.m
% Parameters of the left-ventricle waveform generator based on
% Ferreira et al. (2005), reduced model with constant atrial pressure.
%
% Source:
%   Ferreira A, Chen S, Simaan MA, Boston JR, Antaki JF. A nonlinear
%   state-space model of a combined cardiovascular system and a rotary pump.
%   Proc. 44th IEEE CDC-ECC 2005, Seville, pp. 897-902.
%   Eqs. (1)-(3), (10)-(12), Table II.
%
% Reduced model (z-state form), 3 states:
%   z1 = LVV - V0  [mL]     LVP = E(t)*z1
%   AP             [mmHg]   arterial pressure (compliance C3, drains via R1)
%   QA             [mL/s]   aortic flow (inertance L)
% The atrial pressure is constant (LAP), so the atrial compliance C2
% of the original model is not needed.
% Note: in eq. (12) the element (1,2) is printed as E(t)/R2; the correct
% value is 1/R2 (see eq. 27 and the transformation in eq. 6).
%
% Units: pressure [mmHg], flow [mL/s], volume [mL]
%
% Usage: run AFTER parameters_stergiopulos.m (R1 and LAP are taken from
% there), e.g. in the model's InitFcn callback:
%   parameters_stergiopulos; ferreira_params;

%% Check that the Stergiopulos parameters are loaded
if ~exist('R_periph', 'var') || ~exist('P_venous', 'var') || ~exist('C_total', 'var')
    error(['ferreira_params: run parameters_stergiopulos.m first ' ...
           '(R_periph and P_venous are needed).']);
end

%% ===================== Heart rate and elastance (eqs. 2-3) ======================
HR   = 75;               % [bpm] heart rate
Emax = 2;              % [mmHg/mL] end-systolic elastance (slope of ESPVR)
Emin = 0.06;             % [mmHg/mL] end-diastolic elastance (slope of EDPVR)
                         %   with LAP = 5 mmHg this gives a low filling volume;
                         %   lower it (e.g. 0.04) if pressures come out too low
tc   = 60/HR;            % [s] cardiac period
Tmax = 0.2 + 0.15*tc;    % [s] time to peak elastance (normalises the double-hill function)

%% ===================== Valves and aortic branch (Table II) ======================
R2 = 0.005;              % [mmHg*s/mL]   mitral valve resistance   -> filling phase
R3 = 0.001;              % [mmHg*s/mL]   aortic valve resistance   -> ejection phase
R4 = 0.0398;             % [mmHg*s/mL]   characteristic resistance -> ejection phase
L  = 0.0005;             % [mmHg*s^2/mL] aortic inertance          -> ejection phase

%% ===================== Arterial compartment (generator only) ======================
% Internal arterial load of the generator: it only shapes the ejection.
% The real load is the physical circuit downstream of the chamber.
R1 = R_periph;           % [mmHg*s/mL] systemic resistance = total peripheral
                         %   resistance of Stergiopulos (terminals in parallel)
C3 = C_total;            % [mL/mmHg] total arterial compliance of Stergiopulos
                         %   (sum of segment and terminal compliances);
                         %   paper value (Table II): 1.33
                       
%% ===================== Boundary conditions ======================
LAP = P_venous;          % [mmHg] constant atrial pressure = venous pressure
                         %   of the Stergiopulos model (same downstream node)
V0  = 10;                % [mL] volume at zero pressure (Section IV of the paper)
                         %   only used for the output LVV = z1 + V0;
                         %   TODO: replace with height-sensor calibration

%% ===================== Initial conditions of the integrators ======================
% Only affect the first few beats (transient); steady state is set by
% LAP, Emin/Emax and R1.
%z1_0 = LAP/Emin;         % [mL] end-diastolic filling: LVP = LAP at t = 0
%AP_0 = 80;               % [mmHg] arterial pressure (around diastolic value)
%QA_0 = 0;                % [mL/s] aortic valve closed at t = 0

% Steady-state values at the start of a cycle (end-diastole, aortic valve
% closed), obtained by simulating 10 cycles and reading the final states.
% Valid for: HR 75, Emax 2.0, Emin 0.06, LAP 8, R1 = R_periph, C3 = C_total.
% If any of these parameters changes, re-run and update the values.
z1_0 = 131.9786;            % [mL]   LVV - V0 at end-diastole
AP_0 = 79.3379;             % [mmHg] arterial pressure at end-diastole
QA_0 = 0;                % [mL/s] aortic valve closed

fprintf('Ferreira parameters loaded: HR %d bpm, R1 = %.3f mmHg*s/mL, LAP = %.1f mmHg\n', ...
        HR, R1, LAP);