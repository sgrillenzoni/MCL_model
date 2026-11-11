%% run_ferreira_verification.m
% Standalone verification of the Ferreira et al. (2005) generator.
%   1) runs the nominal case and prints the metrics vs the paper
%   2) (optional) calibrates the constant left atrial pressure LAP, in one of
%      two modes:
%        'SV'  : match the stroke volume of the paper exactly (SV = 67.5 mL)
%        'LSQ' : minimise the sum of squared RELATIVE errors on SBP, DBP, MAP
%                and SV (equal weights). CO is excluded because at fixed HR it
%                is proportional to SV (it would be counted twice), PP because
%                it is SBP - DBP.
%
% Why LAP may be calibrated but Emax may not: Emax is a parameter of the
% paper; the constant LAP is our simplification (the paper has a compliant
% atrium), so the paper gives no value for it. With one parameter only one
% target can be matched: SV is chosen because it is the quantity directly
% set by ventricular filling, and it drives the mean flows of the twin.
% MAP then follows the cardiac output (MAP ~ LAP + CO*R1).
%
% Before running: name the wires 'AoP' (or 'AP') and 'LVV' and enable
% signal logging;
% check the three settings below.

%% ---------------- Settings (CHECK THESE) ----------------
MODEL      = 'ferreira_charlton';
CALIB_MODE = 'SV';
SV_TARGET  = 66.3;            % pressure compared with the paper:
LAP_VAR  = 'LAP';                 % name of the workspace variable used for LAP
P_SIG    = 'AoP';                 % pressure compared with the targets
DO_CALIB   = true;                % true = also calibrate LAP
CALIB_MODE = 'SV';                % 'SV' (match SV) or 'LSQ' (least squares)

SV_TARGET = 66.3;                 % [mL] PWDB baseline 25 y (Charlton 2019)
% Paper values for the least-squares mode: [SBP DBP MAP SV]
REF = [115 75 98 67.5];           % [mmHg mmHg mmHg mL]
W   = [1 1 1 1];              % weights (equal)
STOP_TIME = '20';             % [s] ~25 beats at 75 bpm, enough to reach steady state

%% ---------------- Nominal run ----------------
% Parameters are loaded by the model callbacks (PreLoadFcn / InitFcn).
load_system(MODEL);
HR_run  = evalin('base', 'HR');
LAP_nom = evalin('base', LAP_VAR);

in  = Simulink.SimulationInput(MODEL);
in  = in.setModelParameter('StopTime', STOP_TIME);
out = sim(in);

fprintf('\n=== Nominal case: %s = %.2f mmHg ===\n', LAP_VAR, LAP_nom);
m_nom = ferreira_metrics(out, HR_run, true, P_SIG);

%% ---------------- LAP calibration on SV ----------------
if DO_CALIB
    switch CALIB_MODE
        case 'SV'
            % Residual SV(LAP) - SV_TARGET. SV increases monotonically with
            % LAP (higher filling pressure -> larger EDV): bracketed root.
            svErr   = @(lap) svOf(in, LAP_VAR, lap, HR_run, P_SIG) - SV_TARGET;
            opts    = optimset('TolX', 0.01, 'Display', 'iter');
            LAP_cal = fzero(svErr, [LAP_nom - 3, LAP_nom + 1], opts);
        case 'LSQ'
            % Cost J(LAP) = sum W.*((x - REF)./REF).^2, minimised on a bounded
            % interval (one parameter -> fminbnd, no toolbox needed).
            J       = @(lap) costOf(in, LAP_VAR, lap, HR_run, REF, W, P_SIG);
            opts    = optimset('TolX', 0.01, 'Display', 'iter');
            [LAP_cal, Jmin] = fminbnd(J, LAP_nom - 3, LAP_nom + 1, opts);
            fprintf('Minimum cost J = %.4g (RMS relative error %.1f%%)\n', ...
                    Jmin, 100*sqrt(Jmin/sum(W)));
        otherwise
            error('CALIB_MODE must be ''SV'' or ''LSQ''.');
    end

    out_cal = sim(in.setVariable(LAP_VAR, LAP_cal));
    fprintf('\n=== Calibrated case (%s): %s = %.2f mmHg ===\n', CALIB_MODE, LAP_VAR, LAP_cal);
    m_cal = ferreira_metrics(out_cal, HR_run, true, P_SIG);

    fprintf(['\nNote: setVariable only affects these runs. To use the calibrated ' ...
             'value, write %s = %.2f in the parameter script.\n'], LAP_VAR, LAP_cal);
end

%% ---------------- Local function ----------------
function sv = svOf(in, lapVar, lap, HR, pSig)
% Stroke volume [mL] for a given LAP (one simulation).
out = sim(in.setVariable(lapVar, lap));
m   = ferreira_metrics(out, HR, false, pSig);
sv  = m.SV;
end

function J = costOf(in, lapVar, lap, HR, ref, w, pSig)
% Weighted sum of squared relative errors on [SBP DBP MAP SV].
out = sim(in.setVariable(lapVar, lap));
m   = ferreira_metrics(out, HR, false, pSig);
x   = [m.SBP m.DBP m.MAP m.SV];
J   = sum(w .* ((x - ref)./ref).^2);
end
