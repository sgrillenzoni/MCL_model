%% make_charlton_inflow.m
% Aortic inflow waveform of the PWDB baseline 25-year-old subject, as a
% one-beat lookup table for the twin (flow source at the aortic root).
%
% Uses AorticFlowWave_orig.m from the PWDB code (Charlton et al. 2019, MIT
% licence), with the same inputs as calculate_pwdb_input_parameters.m:
% template shape of Mynard & Smolich 2015, scaled to HR, SV, time of peak
% flow, ejection time and reverse flow volume of the subject.
%
% Requires: parameters_charlton.m run first (CH_HR, CH_SV, CH_PFT, CH_LVET, CH_RFV).
% Output in the workspace:
%   CH_T          [s]     beat period (60/CH_HR)
%   CH_INFLOW_t   [s]     time within the beat, 0..CH_T
%   CH_INFLOW_Q   [mL/s]  aortic flow
%
% In Simulink:  Clock -> mod(u, CH_T) -> 1-D Lookup Table (breakpoints
% CH_INFLOW_t, table CH_INFLOW_Q) -> Simulink-PS Converter (unit A)
% -> Controlled Current Source at the aortic root.

p = struct();
p.wave_type     = 'Mynard';
p.HR            = CH_HR;            % [bpm]
p.SV            = CH_SV;            % [mL]
p.T_Peak_Flow   = 1000*CH_PFT;      % [ms]
p.LVET          = 1000*CH_LVET;     % [ms]
p.Reg_Vol       = CH_RFV;           % [mL]
p.plot_multiple = true;             % do not close open figures
p.do_plot       = false;
p.save_plot     = false;
p.plot_name     = 'AorticFlowWave_plot';
p.file_path     = tempdir;     % avoids a folder dialog in AorticFlowWave_orig

w = AorticFlowWave_orig(p);

CH_T        = 60/CH_HR;
CH_INFLOW_t = w.t(:)';
CH_INFLOW_Q = 1e6*w.v(:)';          % m^3/s -> mL/s
% close the beat exactly at CH_T, so that the table is periodic
if CH_INFLOW_t(end) < CH_T - 1e-9
    CH_INFLOW_t(end+1) = CH_T;
    CH_INFLOW_Q(end+1) = CH_INFLOW_Q(1);
else
    CH_INFLOW_t(end) = CH_T;
end

fprintf('PWDB inflow: HR %.1f bpm, SV %.1f mL (target %.1f), peak %.0f mL/s at %.0f ms, CO %.2f L/min\n', ...
        CH_HR, trapz(CH_INFLOW_t, CH_INFLOW_Q), CH_SV, max(CH_INFLOW_Q), ...
        1000*CH_INFLOW_t(find(CH_INFLOW_Q == max(CH_INFLOW_Q), 1)), ...
        trapz(CH_INFLOW_t, CH_INFLOW_Q)*CH_HR/1000);
clear p w
