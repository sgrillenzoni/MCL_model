function m = ferreira_metrics(out, HR, doPrint, pSig)
%FERREIRA_METRICS  Haemodynamic metrics on the last beat of a simulation.
%   m = ferreira_metrics(out, HR)                 computes and prints (on AoP)
%   m = ferreira_metrics(out, HR, false)          computes without printing
%   m = ferreira_metrics(out, HR, true, 'AP')     uses another pressure signal
%
%   out : Simulink.SimulationOutput returned by sim() (signal logging on)
%   HR  : heart rate [bpm]
%   pSig: name of the logged pressure signal (default 'AoP').
%         'AoP' = aortic root pressure (after the valve, before R4 and L): the
%                 quantity reported in Ferreira 2005, Fig. 4 -> use it to
%                 compare with the paper.
%         'AP'  = arterial node across C3 (model state).
%   m   : struct with fields SBP, DBP, MAP, PP [mmHg], SV [mL], CO [L/min]
%
%   Requires the logged signals to be NAMED (signal name on the wire, not the
%   Outport name): the pressure signal pSig and 'LVV'.

if nargin < 3 || isempty(doPrint), doPrint = true;  end
if nargin < 4 || isempty(pSig),    pSig    = 'AoP'; end

ls  = out.logsout;
AP  = getSignal(ls, pSig);      % pressure used for SBP/DBP/MAP
LVV = getSignal(ls, 'LVV');

T  = 60/HR;                       % beat period [s]
t1 = AP.Time(end) - T;            % start of the last beat
iA = AP.Time  >= t1;
iV = LVV.Time >= t1;

m.SBP = max(AP.Data(iA));
m.DBP = min(AP.Data(iA));
m.MAP = trapz(AP.Time(iA), AP.Data(iA)) / (AP.Time(end) - t1);  % time-weighted
m.PP  = m.SBP - m.DBP;
m.SV  = max(LVV.Data(iV)) - min(LVV.Data(iV));
m.CO  = m.SV * HR / 1000;

if doPrint
    % Reference: Ferreira et al. 2005, nominal case (Fig. 4), HR = 75 bpm
    fprintf('Pressure signal: %s\n', pSig);
    fprintf('SBP/DBP %6.1f / %5.1f mmHg   (paper 115 / 75)\n', m.SBP, m.DBP);
    fprintf('MAP %6.1f mmHg           (paper 98)    %+5.1f%%\n', m.MAP, 100*(m.MAP/98   - 1));
    fprintf('PP  %6.1f mmHg           (paper 40)    %+5.1f%%\n', m.PP,  100*(m.PP/40    - 1));
    fprintf('SV  %6.1f mL             (paper 67.5)  %+5.1f%%\n', m.SV,  100*(m.SV/67.5  - 1));
    fprintf('CO  %6.2f L/min          (paper 5.06)  %+5.1f%%\n', m.CO,  100*(m.CO/5.06  - 1));
end
end

function v = getSignal(ls, name)
% Find a logged signal by name; if missing, list the available names.
e = ls.getElement(name);
if isempty(e)
    names = cell(1, ls.numElements);
    for k = 1:ls.numElements
        names{k} = ls.getElement(k).Name;
    end
    error(['ferreira_metrics: logged signal "%s" not found. Available: %s\n' ...
           'Name the wire (double-click on it) and enable logging.'], ...
           name, strjoin(strcat('"', names, '"'), ', '));
end
if isa(e, 'Simulink.SimulationData.Dataset')   % several signals with same name
    e = e.getElement(1);
end
v = e.Values;
end
