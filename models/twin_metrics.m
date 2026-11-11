function m = twin_metrics(out, HR, names)
%TWIN_METRICS  Global haemodynamics and flow split of the digital twin,
%   computed on the last beat of a simulation.
%
%   m = twin_metrics(out, HR)          uses the default signal names below
%   m = twin_metrics(out, HR, names)   names = struct with fields
%                                      P, Q_LV, Q_upper, Q_abdominal, Q_lower
%                                      and optionally P_sites (cell array)
%
%   out : Simulink.SimulationOutput returned by sim() (signal logging on)
%   HR  : heart rate [bpm]
%
%   Logged signals needed (named on the wire):
%     P            aortic pressure at the outlet of the generator [mmHg]
%     Q_ao         flow ejected by the ventricle [mL/s]
%     Q_upper      mean flow into the upper body chamber [mL/s]
%     Q_abdominal  mean flow into the abdominal chamber [mL/s]
%     Q_lower      mean flow into the lower body chamber [mL/s]
%     P_sites      (optional) pressures to compare site by site, first one
%                  = reference for the delay. Default: AoP and the pressure
%                  targets of the three chambers (the same quantities that
%                  are measured on the bench).
%
%   Reference values printed for comparison:
%   - expected split of the Stergiopulos network at zero frequency, from its
%     resistances (terminal + own segment, in parallel): checks the
%     implementation of the network (verification)
%   - resting distribution of cardiac output, Williams & Leggett (1989):
%     kidneys 19%, liver total 25% (hepatic artery + gut, spleen, pancreas),
%     brain 12% (validation; not a one-to-one match with the chambers)

if nargin < 3 || isempty(names)
    names = struct('P', 'AoP', 'Q_LV', 'Q_LV', 'Q_upper', 'Q_upper', ...
                   'Q_abdominal', 'Q_abdominal', 'Q_lower', 'Q_lower');
end

if ~isfield(names, 'P_sites') || isempty(names.P_sites)
    names.P_sites = {names.P, 'P_target_upperbody', 'P_target_abdominal', 'P_target_lowerbody'};
end

ls = out.logsout;
T  = 60/HR;

% --- Pressure ---
P  = getSignal(ls, names.P);
[tP, yP] = lastBeat(P, T);
m.SBP = max(yP);
m.DBP = min(yP);
m.MAP = beatMean(tP, yP);
m.PP  = m.SBP - m.DBP;

% --- Flows: mean over the last beat ---
m.CO_mLs = beatMean2(getSignal(ls, names.Q_Ao), T);          % [mL/s]
m.CO     = m.CO_mLs * 60/1000;                                 % [L/min]
m.SV     = m.CO_mLs * 60/HR;                                   % [mL]
q = [beatMean2(getSignal(ls, names.Q_upper),     T), ...
     beatMean2(getSignal(ls, names.Q_abdominal), T), ...
     beatMean2(getSignal(ls, names.Q_lower),     T)];
m.Q_branch = q;                                                % [mL/s]
m.split    = 100 * q / sum(q);                                 % [% of total]
m.balance  = 100 * (sum(q) / m.CO_mLs - 1);                    % mass check [%]

% --- Print ---
fprintf('\n--- Global haemodynamics (last beat) ---\n');
fprintf('Aortic pressure %6.1f / %5.1f mmHg   (normal adult at rest ~120/80)\n', m.SBP, m.DBP);
fprintf('MAP             %6.1f mmHg\n', m.MAP);
fprintf('Pulse pressure  %6.1f mmHg\n', m.PP);
fprintf('SV              %6.1f mL\n', m.SV);
fprintf('CO              %6.2f L/min          (normal adult at rest ~5)\n', m.CO);
fprintf('Mass balance    %+6.2f %%  (sum of branch flows vs ejected flow; should be ~0)\n', m.balance);

ref_net = [18.5 62.7 18.7];   % Stergiopulos network, zero-frequency split [%]
lbl = {'Upper body', 'Abdominal', 'Lower body'};
fprintf('\n--- Flow split (mean over last beat) ---\n');
fprintf('%-12s %10s %10s %16s\n', 'Chamber', 'Q [mL/s]', 'Share [%]', 'Expected net [%]');
for k = 1:3
    fprintf('%-12s %10.1f %10.1f %16.1f\n', lbl{k}, q(k), m.split(k), ref_net(k));
end
fprintf(['Physiology (Williams & Leggett 1989): kidneys 19%%, liver total 25%%, ' ...
         'brain 12%% of cardiac output.\n']);

% --- Pressures site by site ---
% Delay = arrival of the pressure wave at each site minus at the reference
% (first in the list). Arrival = first crossing of 20% of the pulse amplitude
% on the systolic upstroke (see upstrokeTime). Two beats are used so that the
% upstroke is never cut by the window edge.
FOOT_FRAC = 0.2;                          % fraction of the pulse amplitude
sites = names.P_sites;
nS    = numel(sites);
m.P_sites = struct('name', sites, 'SBP', [], 'DBP', [], 'MAP', [], 'PP', [], 'delay_ms', []);
tUp   = zeros(1, nS);
for k = 1:nS
    ts = getSignal(ls, sites{k});
    [t, y] = lastBeat(ts, T);
    m.P_sites(k).SBP = max(y);
    m.P_sites(k).DBP = min(y);
    m.P_sites(k).MAP = beatMean(t, y);
    m.P_sites(k).PP  = m.P_sites(k).SBP - m.P_sites(k).DBP;
    tUp(k) = upstrokeTime(ts, T, FOOT_FRAC);
end
for k = 1:nS
    d = mod(tUp(k) - tUp(1), T);
    if d > T/2, d = d - T; end
    m.P_sites(k).delay_ms = 1000 * d;
end

fprintf('\n--- Pressures by site (last beat) ---\n');
fprintf('%-22s %7s %7s %7s %7s %11s\n', 'Site', 'SBP', 'DBP', 'MAP', 'PP', 'Delay [ms]');
for k = 1:nS
    fprintf('%-22s %7.1f %7.1f %7.1f %7.1f %11.0f\n', sites{k}, m.P_sites(k).SBP, ...
            m.P_sites(k).DBP, m.P_sites(k).MAP, m.P_sites(k).PP, m.P_sites(k).delay_ms);
end
fprintf('(pressures in mmHg; delay of the wave foot (20%% of pulse amplitude) vs %s)\n', sites{1});
end

% ---------------- helpers ----------------
function [t, y] = lastBeat(ts, T)
t1 = ts.Time(end) - T;
i  = ts.Time >= t1;
t  = ts.Time(i);
y  = ts.Data(i);
end

function v = beatMean(t, y)
v = trapz(t, y) / (t(end) - t(1));     % time-weighted (works with variable step)
end

function v = beatMean2(ts, T)
[t, y] = lastBeat(ts, T);
v = beatMean(t, y);
end

function tu = upstrokeTime(ts, T, frac)
% Arrival time of the pressure wave: first crossing, on the systolic upstroke,
% of the level  P_min + frac*(P_max - P_min)  (frac = 0.2), searched after the
% diastolic minimum of the second-to-last beat, with linear interpolation.
% A fixed-fraction threshold is more robust than max dP/dt: it does not pick
% the step of AoP at valve closure, nor an early point on a flat diastolic
% minimum, and it is close to the "foot" used for pulse wave velocity.
tEnd = ts.Time(end);
i1   = ts.Time >= tEnd - 2*T & ts.Time < tEnd - T;     % second-to-last beat
t1   = ts.Time(i1);  y1 = ts.Data(i1);
[yMin, iMin] = min(y1);
yMax = max(y1);
lev  = yMin + frac*(yMax - yMin);
i2   = find(ts.Time >= t1(iMin) & ts.Time <= t1(iMin) + T);
t2   = ts.Time(i2);  y2 = ts.Data(i2);
j    = find(y2(1:end-1) < lev & y2(2:end) >= lev, 1, 'first');
if isempty(j)
    tu = NaN;
else
    tu = t2(j) + (lev - y2(j)) * (t2(j+1) - t2(j)) / (y2(j+1) - y2(j));
end
end

function v = getSignal(ls, name)
e = ls.getElement(name);
if isempty(e)
    n = cell(1, ls.numElements);
    for k = 1:ls.numElements, n{k} = ls.getElement(k).Name; end
    error(['twin_metrics: logged signal "%s" not found. Available: %s\n' ...
           'Name the wire (double-click on it) and enable logging.'], ...
           name, strjoin(strcat('"', n, '"'), ', '));
end
if isa(e, 'Simulink.SimulationData.Dataset'), e = e.getElement(1); end
v = e.Values;
end
