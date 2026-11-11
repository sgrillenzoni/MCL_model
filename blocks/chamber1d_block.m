function P_target = chamber1d_block(Q_meas, S)
%CHAMBER1D_BLOCK  Body of the MATLAB Function block of one Windkessel chamber.
%   P_target = chamber1d_block(Q_meas, S)
%     Q_meas   [mL/s]  flow measured into the chamber
%     S                sub-tree from charlton1d_setup (S_upper / S_abdominal /
%                      S_lower), passed as a NON-TUNABLE parameter
%     P_target [mmHg]  pressure the chamber must reproduce
%   Block sample time = S.dt (1e-3 s). One solver step per call.
%#codegen
persistent A Q Pn Pc
if isempty(A)                                   % rest at P = DBP
    A  = S.A0;  Q = zeros(S.nF,1);
    Pn = S.Pref*ones(S.nNodes,1);  Pc = S.Pref*ones(numel(S.tSeg),1);
end
[A, Q, Pn, Pc] = charlton1d_step(S, A, Q, Pn, Pc, 1e-6*Q_meas);
P_target = Pn(S.root)/133.322;
end
