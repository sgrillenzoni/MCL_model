function P_target = twin_chambers_block(Q_meas, S_upper, S_abdominal, S_lower)
%TWIN_CHAMBERS_BLOCK  Body of the MATLAB Function block "Twin1D".
%   One block for the three Windkessel chambers.
%   Q_meas   [mL/s] 3x1 flow measured INTO each chamber (positive entering):
%                   [upper; abdominal; lower]
%   S_upper, S_abdominal, S_lower: 1D sub-trees from parameters_chamber1d.m,
%                   NON-TUNABLE parameters of the block
%   P_target [mmHg] 3x1 pressure each chamber must reproduce
%   Block sample time 1e-3 s (= S_*.dt), one solver step per call.
%#codegen
persistent Au Qu Pnu Pcu  Aa Qa Pna Pca  Al Ql Pnl Pcl
if isempty(Au)                                  % rest at P = DBP
    Au = S_upper.A0;      Qu = zeros(S_upper.nF,1);
    Pnu = S_upper.Pref*ones(S_upper.nNodes,1);      Pcu = S_upper.Pref*ones(numel(S_upper.tSeg),1);
    Aa = S_abdominal.A0;  Qa = zeros(S_abdominal.nF,1);
    Pna = S_abdominal.Pref*ones(S_abdominal.nNodes,1); Pca = S_abdominal.Pref*ones(numel(S_abdominal.tSeg),1);
    Al = S_lower.A0;      Ql = zeros(S_lower.nF,1);
    Pnl = S_lower.Pref*ones(S_lower.nNodes,1);      Pcl = S_lower.Pref*ones(numel(S_lower.tSeg),1);
end
[Au, Qu, Pnu, Pcu] = charlton1d_step(S_upper,     Au, Qu, Pnu, Pcu, 1e-6*Q_meas(1));
[Aa, Qa, Pna, Pca] = charlton1d_step(S_abdominal, Aa, Qa, Pna, Pca, 1e-6*Q_meas(2));
[Al, Ql, Pnl, Pcl] = charlton1d_step(S_lower,     Al, Ql, Pnl, Pcl, 1e-6*Q_meas(3));
P_target = [Pnu(S_upper.root); Pna(S_abdominal.root); Pnl(S_lower.root)]/133.322;
end
