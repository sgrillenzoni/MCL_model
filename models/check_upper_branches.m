%% check_upper_branches.m
% Mean, min and max flow into each branch of the virtual upper-body tree,
% read from the Simscape log, plus the share of cardiac output.
%
% BEFORE RUNNING (once):
%   Model Settings > Simscape > Data Logging > Log simulation data = All
%   (workspace variable name: simlog)
% THEN:
%   1) run the simulation (Run button, or out = sim('model_name'))
%   2) run this script
% Needs in the workspace: out (simulation output) and HR.

%% ---------------- Settings ----------------
% Path of the tree inside the log: block names, spaces become '_'.
% If it fails, explore with  out.simlog.childIds  and go down one level at a time.
TREE        = 'Upperbody_Windkessel1.UpperBody_tree';
NBEATS_PLOT = 2;          % beats shown in the plot

%% ---------------- Get the tree ----------------
t = out.simlog;
for p = strsplit(TREE, '.')
    t = t.(p{1});
end
T = 60/HR;                % beat period [s]

%% ---------------- Branch flows (first cell of each branch) ----------------
% The 5 roots leave the upper-body node: their sum is Q_upper.
% vertebral_right branches off subclavianA_right: used to get the right arm.
roots = {'commonCarotid_right', 'commonCarotid_left', 'vertebral_left', ...
         'subclavianA_right',   'subclavianB_left'};
extra = {'vertebral_right'};
names = [roots extra];

S = struct();
for k = 1:numel(names)
    s  = t.(['R_' names{k} '_1']).i.series;   % current = flow [mL/s]
    tt = s.time;  q = s.values;
    i  = tt >= tt(end) - T;                   % last beat
    S.(names{k}).mean = trapz(tt(i), q(i)) / (tt(end) - tt(find(i, 1)));
    S.(names{k}).min  = min(q(i));
    S.(names{k}).max  = max(q(i));
    S.(names{k}).t    = tt;
    S.(names{k}).q    = q;
end

%% ---------------- Cardiac output (logged signal Q_Ao) ----------------
CO = NaN;
try
    e = out.logsout.getElement('Q_Ao');
    if isa(e, 'Simulink.SimulationData.Dataset'), e = e.getElement(1); end
    tq = e.Values.Time;  qa = e.Values.Data;
    i  = tq >= tq(end) - T;
    CO = trapz(tq(i), qa(i)) / (tq(end) - tq(find(i, 1)));   % [mL/s]
catch
    warning('Q_Ao not found in logsout: shares of CO not computed.');
end

%% ---------------- Table ----------------
fprintf('\n%-22s %8s %8s %8s %10s\n', 'Branch', 'mean', 'min', 'max', 'back [%]');
for k = 1:numel(names)
    b = S.(names{k});
    bf = max(0, -b.min) / b.max * 100;        % backflow as % of peak
    fprintf('%-22s %8.2f %8.2f %8.2f %10.0f\n', names{k}, b.mean, b.min, b.max, bf);
end

armR  = S.subclavianA_right.mean - S.vertebral_right.mean;
armL  = S.subclavianB_left.mean;
upper = 0;
for k = 1:numel(roots), upper = upper + S.(roots{k}).mean; end

fprintf('\n--- Derived (mean flows, mL/s) ---\n');
fprintf('Right arm (subclavian A R - vertebral R) %6.2f\n', armR);
fprintf('Left arm  (subclavian B L)               %6.2f\n', armL);
fprintf('Vertebrals (R + L)                       %6.2f\n', S.vertebral_right.mean + S.vertebral_left.mean);
fprintf('Carotids  (R + L)                        %6.2f\n', S.commonCarotid_right.mean + S.commonCarotid_left.mean);
fprintf('Upper total (sum of the 5 roots)         %6.2f\n', upper);
if ~isnan(CO)
    fprintf('Cardiac output                           %6.2f  (%.2f L/min)\n', CO, CO*60/1000);
    fprintf('Upper share of CO                        %6.1f %%\n', 100*upper/CO);
end

%% ---------------- Plot: last beats of the 5 roots ----------------
figure; hold on; grid on
for k = 1:numel(roots)
    b = S.(roots{k});
    i = b.t >= b.t(end) - NBEATS_PLOT*T;
    plot(b.t(i), b.q(i), 'DisplayName', roots{k}, 'LineWidth', 1.2);
end
yline(0, 'k:', 'HandleVisibility', 'off');
legend('Interpreter', 'none', 'Location', 'best');
xlabel('t [s]'); ylabel('Q [mL/s]');
title('Upper body: flow into each branch');
