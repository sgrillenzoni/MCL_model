%% parameters_backflowcontrol.m
% Parameters for the return-flow (backflow) control of the Windkessel
% chambers: the pumps that bring the water back from each chamber to the
% heart chamber.
%
% Control structure (per chamber):
%   [1] CycleMean    one-beat moving average of the sensor signals
%   [2] FeedForward  mean inlet flow -> pump command (pump LUT)
%   [3] Feedback     slow PI on the level error (h_mean - h_ref)
%   [4] Calibration  captures h_ref at steady state, heart beating
%   [5] Safety       min/max level, saturation, rate limit
%
% Units:
%   time [s], heart rate [bpm], flow [mL/s], level [m], pump command [V]
%   (level in m, as in the lab height-control model, h_hp [m])
%
% Usage: run this script before the simulation, or call it from the model's
% PreLoadFcn / InitFcn callback (Model Properties > Callbacks).
% Blocks use these variable names directly (no masks).
%
% Currently set up for ONE chamber. When adding the other chambers, give the
% chamber-specific parameters a suffix (e.g. LUT_Q_upper, LUT_Q_abdominal,
% LUT_Q_lower); the timing parameters (HR, Ts, N) stay shared.

%% ===================== Timing ======================

% Heart rate [bpm]: NOT defined here. It is taken from the Ferreira parameter
% script, so that the one-beat average always uses the same heart rate as the
% pressure generator (a mismatch leaves a residual ripple).
% -> The Ferreira script must be run BEFORE this one (same order in the
%    model callbacks).
if ~exist('HR', 'var')
    error(['parameters_backflowcontrol: HR not defined. ' ...
           'Run the Ferreira parameter script before this one.']);
end

% Base sample time of the real-time model (dSPACE solver, fixed-step) [s].
% Sensor signals and the one-beat average run at this rate.
% The whole return control runs at this single rate (no Rate Transition):
% times in the later blocks are expressed in seconds, not in samples.
Ts = 0.0002;

%% ============ [1] CycleMean - one-beat moving average ============

% Number of samples in one cardiac period (window length of the average).
N = round(60/HR/Ts);

% The window covers exactly one beat only if 60/HR/Ts is an integer;
% otherwise a small ripple remains in the averaged signal.
if abs(60/HR/Ts - N) > 1e-6
    warning(['parameters_backflowcontrol: 60/HR/Ts is not an integer, ' ...
             'the one-beat window is not exact (small residual ripple).']);
end

%% ============ [2] FeedForward - pump look-up table ============

% Pump characteristic in the real circuit configuration:
% for each flow LUT_Q, the command LUT_u that produces it.
% Used as an inverse LUT: mean inlet flow -> pump command.
%   - LUT_Q must be strictly increasing (pump curve monotonic in the used range)
%   - first point (0, 0); second point (Q_min, u_min) = dead zone of the pump,
%     i.e. the minimum command at which the pump actually delivers flow
%   - outside the measured range the LUT block clips (no extrapolation)
%
% Which LUT to use:
%   'mock'      placeholder values below, for testing the blocks
%   'measured'  measured curve of the Jabsco 23680-4103 in the real
%               configuration, loaded from LUT_FILE (variables LUT_Q, LUT_u)
LUT_SOURCE = 'mock';
LUT_FILE   = 'LUT_pump_measured.mat';

switch LUT_SOURCE
    case 'mock'
        % !!! PLACEHOLDER VALUES, for testing only - NOT the real pump.
        LUT_Q = [0  5    10   20   40];    % [mL/s]
        LUT_u = [0  1.5  2.5  4.5  8 ];    % [V]
        warning('parameters_backflowcontrol: using MOCK pump LUT (test values).');
    case 'measured'
        if ~exist(LUT_FILE, 'file')
            error('parameters_backflowcontrol: measured LUT file "%s" not found.', LUT_FILE);
        end
        S = load(LUT_FILE, 'LUT_Q', 'LUT_u');
        LUT_Q = S.LUT_Q(:)';    % row vectors
        LUT_u = S.LUT_u(:)';
        clear S
    otherwise
        error('parameters_backflowcontrol: LUT_SOURCE must be ''mock'' or ''measured''.');
end

if numel(LUT_Q) ~= numel(LUT_u)
    error('parameters_backflowcontrol: LUT_Q and LUT_u must have the same length.');
end
if any(diff(LUT_Q) <= 0)
    error('parameters_backflowcontrol: LUT_Q must be strictly increasing.');
end

%% ============ [3] Feedback - slow PI on the chamber level ============

% Chamber geometry (same tanks as the twin circuit).
tank_radius = 0.06;                     % [m]
tank_height = 0.14;                     % [m]
A_tank = pi*tank_radius^2 * 1e6;        % [mL/m]  water volume per metre of level
% (The floater displaces a constant volume, so it does not change A_tank.)

% Mean pump gain k_pump [mL/s per V]: slope of the LUT between the end of the
% dead zone (2nd point) and the last point. A rough mean value is enough:
% an error here only makes the feedback somewhat faster or slower, and the
% integral term still brings the level exactly to h_ref.
k_pump = (LUT_Q(end) - LUT_Q(2)) / (LUT_u(end) - LUT_u(2));

% Desired time constant of the level feedback [s]: tens of seconds, i.e. much
% slower than one beat (the one-beat mean has a delay of half a beat) and
% slow enough not to fight the feedforward, which does most of the work.
% TUNING ON THE BENCH: 30 s is a conservative starting value. In simulation
% (pump 10% weaker than the LUT) the level peaks at +4 mm and returns within
% 0.5 mm of h_ref in ~270 s. Faster options tested in simulation:
%   tau_fb = 15 s -> peak +2.2 mm, within 0.5 mm in ~110 s
%   tau_fb = 10 s -> peak +1.6 mm, within 0.5 mm in ~60 s
% Decrease step by step on the bench: a faster loop also reacts more to
% sensor noise and surface waves. Change only tau_fb (Ti = 4*tau_fb);
% faster integral alone (Ti = 2*tau_fb) gives a small undershoot.
tau_fb = 30;

% PI gains, parallel form  u_fb = Kp_fb*e + Ki_fb*integral(e),  e = h_mean - h_ref
%   Kp_fb = A/(k*tau): the chamber is itself an integrator of the flow
%                      imbalance, so this gives a first-order response with
%                      time constant tau_fb.
%   Ti_fb = 4*tau_fb:  integral slower than the proportional action, otherwise
%                      the two integrators (chamber + PI) make the level oscillate.
Kp_fb = A_tank / (k_pump * tau_fb);     % [V/m]
Ti_fb = 4 * tau_fb;                     % [s]
Ki_fb = Kp_fb / Ti_fb;                  % [V/(m*s)]
% For reference, the lab height-control PI (no feedforward): P = 120, I = 150.

% Limited authority of the feedback [V]: the correction can move the command
% only +/- U_FB_MAX around the feedforward value. If the level sensor fails,
% the pump stays close to the (sensible) feedforward command. Anti-windup
% (clamping) acts on this limit.
U_FB_MAX = 2;

% Setpoint saturation [m]: h_ref captured by the calibration is limited to
% this range, so an absurd value (sensor error during settling) is never
% tracked. !!! Placeholder limits - check against the real chamber.
H_REF_MIN = 0.03;
H_REF_MAX = 0.12;

%% Summary
fprintf('Backflow control parameters loaded: HR = %g bpm, Ts = %g s, N = %d samples/beat, LUT = %s.\n', ...
        HR, Ts, N, LUT_SOURCE);
fprintf('Level feedback: k_pump = %.2f mL/s/V, Kp_fb = %.1f V/m, Ti_fb = %g s, Ki_fb = %.3f V/(m*s).\n', ...
        k_pump, Kp_fb, Ti_fb, Ki_fb);
