% Loads all parameters needed by the MCL model.
% Add one line here for each new parameter file.
parameters_stergiopulos;
parameters_ferreira;

MODE = 1; % 1 = numerical, 2 = testbench
%% Butterworth low-pass for flow display (2nd order, analog)
fc_flow = 10;                          % [Hz] cutoff frequency
wc      = 2*pi*fc_flow;                % [rad/s]
num_bw  = wc^2;
den_bw  = [1, sqrt(2)*wc, wc^2];