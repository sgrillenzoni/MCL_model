function y = twin_full_block(LVP, S_full)

%TWIN_FULL_BLOCK  Body of the MATLAB Function block "TwinFull1D": the whole
%   PWDB network (aorta of Charlton + upper, abdominal and lower trees).
%   LVP   [mmHg]  left-ventricular pressure (Ferreira generator), upstream of the valve
%   y     8x1:    [P_root; Q_upper; Q_abdominal; Q_lower; P_upper; P_abdominal; P_lower; Q_ao]
%   S_full        full network from parameters_chamber1d.m, NON-TUNABLE parameter
%   Block sample time 1e-3 s (= S_full.dt).
%#codegen
persistent A Q Pn Pc
if isempty(A)                                   % rest at P = DBP
    A  = S_full.A0;  Q = zeros(S_full.nF,1);
    Pn = S_full.Pref*ones(S_full.nNodes,1);  Pc = S_full.Pref*ones(numel(S_full.tSeg),1);
end
[A, Q, Pn, Pc, Qv] = charlton1d_step(S_full, A, Q, Pn, Pc, 0, ...
                                     LVP*133.322, S_full.Rv);
y = zeros(8,1);
y(1) = Pn(S_full.root)/133.322;
y(2) = 1e6*sum(Q(S_full.fUp));
y(3) = 1e6*sum(Q(S_full.fAb));
y(4) = 1e6*sum(Q(S_full.fLo));
y(5) = mean(cellP(S_full, A, S_full.cUp));
y(6) = mean(cellP(S_full, A, S_full.cAb));
y(7) = mean(cellP(S_full, A, S_full.cLo));
y(8) = 1e6*Qv;                                   % [mL/s]
end

function p = cellP(S, A, c)                     % [mmHg] elastic pressure in cells c
p = (S.Pref + S.beta(c)./S.A0(c).*(sqrt(A(c)) - S.sqA0(c)))/133.322;
end
