function [A, Q, Pn, Pc, Qv] = charlton1d_step(S, A, Q, Pn, Pc, Qin, Pin, Rv)
%CHARLTON1D_STEP  One explicit fixed time step of the 1D network
%   (staggered finite volumes, symplectic Euler; see charlton1d_setup.m).
%
%   [A, Q, Pn, Pc] = charlton1d_step(S, A, Q, Pn, Pc, Qin)
%     A   [m^2]   cell areas                    (S.nC x 1)
%     Q   [m^3/s] face flows                    (S.nF x 1)
%     Pn  [Pa]    network-node pressures        (S.nNodes x 1)
%     Pc  [Pa]    Windkessel compliance pressures (numel(S.tSeg) x 1)
%     Qin [m^3/s] inflow into the root node
%   Initial state (rest at P = DBP):
%     A = S.A0, Q = zeros(S.nF,1), Pn = S.Pref*ones(S.nNodes,1),
%     Pc = S.Pref*ones(numel(S.tSeg),1)
%
%   Order: momentum on the faces (old areas) -> node pressures (linear) ->
%   wall viscosity (implicit, tridiagonal) -> mass on the cells (new flows)
%   -> Windkessels (forward Euler).
%   Written with plain loops and no sparse matrices: compatible with code
%   generation (MATLAB Function block, Simulink Coder / dSPACE).

if nargin < 8, Pin = 0; Rv = 0; end             % no valve: flow inlet as before

dt = S.dt; rho = S.rho; Kr = S.Kr;
P = S.Pref + S.beta./S.A0.*(sqrt(A) - S.sqA0);   % elastic pressure in the cells
Qn = Q;

% ---- interior faces: momentum ----
f = S.fInt; L = S.fL(f); R = S.fR(f); dx = S.fdx(f);
Qc = (Q(S.cLF) + Q(S.cRF))/2;                     % flow at cell centres
Fc = S.alpha*Qc.^2./A;                            % convective flux
Af = (A(L) + A(R))/2;
Qn(f) = (Q(f) + dt*( -Af/rho.*(P(R) - P(L))./dx - (Fc(R) - Fc(L))./dx )) ...
        ./ (1 + dt*Kr./Af);

% ---- end faces + network nodes: Q_end = a - s*g*Pn, one linear eq. per node ----
e = S.eF; s = S.eS; c = S.eC;
k = A(c)./(rho*S.fdx(e)/2);
D = 1 + dt*Kr./A(c);
a = (Q(e) + s.*dt.*k.*P(c))./D;
g = dt.*k./D;
num = zeros(S.nNodes,1);  den = zeros(S.nNodes,1);
for j = 1:numel(e)
    num(S.eNode(j)) = num(S.eNode(j)) + s(j)*a(j);
    den(S.eNode(j)) = den(S.eNode(j)) + g(j);
end
num(S.root) = num(S.root) + Qin;
for j = 1:numel(S.tNode)                          % (a node may carry >1 Windkessel)
    num(S.tNode(j)) = num(S.tNode(j)) + Pc(j)/S.R1(j);
    den(S.tNode(j)) = den(S.tNode(j)) + 1/S.R1(j);
end
Pn = num./den;
Qv = 0;
r = S.root;
if Rv > 0 && Pin > Pn(r)                         % aortic valve open (ideal diode + Rv)
    Pn(r) = (num(r) + Pin/Rv)/(den(r) + 1/Rv);
    Qv = (Pin - Pn(r))/Rv;                       % [m^3/s] flow through the valve
end
Qn(e) = a - s.*g.*Pn(S.eNode);

% ---- wall viscosity: implicit diffusion, Thomas algorithm (factorised in setup) ----
nF = S.nF;
y = zeros(nF,1);
y(1) = Qn(1)*S.tm(1);
for i = 2:nF
    y(i) = (Qn(i) - S.tdL(i)*y(i-1))*S.tm(i);
end
for i = nF-1:-1:1
    y(i) = y(i) - S.tc(i)*y(i+1);
end
Q = y;

% ---- cells: mass with the new flows (exactly conservative) ----
A = A - dt./S.dxc.*(Q(S.cRF) - Q(S.cLF));

% ---- Windkessels (forward Euler) ----
Q1 = (Pn(S.tNode) - Pc)./S.R1;
Pc = Pc + dt./S.CT.*(Q1 - (Pc - S.Pout)./S.R2);
end
