function S = charlton1d_setup(geoFile, prm, dt, dxMax, region, stubF)
%CHARLTON1D_SETUP  Precompute the 1D arterial network of the PWDB
%   (Charlton et al. 2019, baseline 25 y), or one chamber sub-tree of it,
%   for the explicit fixed-step solver charlton1d_step.m.
%
%   S = charlton1d_setup(geoFile, prm, dt, dxMax)                 whole network
%   S = charlton1d_setup(geoFile, prm, dt, dxMax, region, stubF)  sub-tree
%     geoFile : 'pwdb_geo0001.csv' (116 segments)
%     prm     : struct with rho [kg/m^3], mu [Pa s], alpha [-], k1 k2 k3
%               (stiffness, PWDB units), b0 b1 (wall viscosity, PWDB units),
%               Pref [mmHg] (pressure at which the CSV radii hold = DBP),
%               Pout [mmHg] (outflow pressure of the Windkessels)
%     dt      : time step [s] (fixed)
%     dxMax   : maximum cell length [m] (e.g. 0.02)
%     region  : 'full' (default) | 'upper' | 'abdominal' | 'lower'
%     stubF   : (upper only) physical fraction of the four phantom outlets
%               [R subclavian A (4), R common carotid (5), L common carotid (15),
%                L subclavian A (19)], e.g. [F_STUB_SUBCLAVIAN_R ...]
%
%   The inflow Qin of charlton1d_step enters the ROOT node:
%     full      : inlet of the ascending aorta (segment 1)
%     upper     : one node joining the virtual remainders of segments 4, 5,
%                 15, 19 (the part after the phantom stub). A remainder that
%                 is fully physical is skipped and its children start at the root.
%     abdominal : one node joining the inlets of 26 (intercostals), 29 (celiac),
%                 34 (sup. mesenteric), 36, 38 (renals), 40 (inf. mesenteric)
%     lower     : one node joining the inlets of 42, 43 (common iliacs)
%   The pressure of the root node, S_out = st.Pn(S.root), is the chamber
%   pressure target.
%
%   Model (as in the PWDB / Nektar1D input files, PWDB code
%   create_input_files_for_a_simulation.m):
%     dA/dt + dQ/dx = 0
%     dQ/dt + d(alpha Q^2/A)/dx + (A/rho) dP/dx = -Kr Q/A,
%     Kr = 2*pi*mu*alpha/(rho*(alpha-1)) = 8*pi*mu/rho for alpha = 4/3
%     P = Pref + beta/A0*(sqrt(A) - sqrt(A0)) + Gamma/(A0*sqrt(A)) dA/dt
%     A0 = pi r^2 at P = Pref (= DBP),  r linear along each segment
%     Eh = 0.1*(k1*exp(k2*100*r) + k3)*r,  beta = 4/3*sqrt(pi)*Eh
%     optional prm.Kstiff (one value per geometry segment, default 1):
%     Eh -> Kstiff*Eh, i.e. local compliance / Kstiff, wave speed * sqrt(Kstiff)
%     Gamma = (b0 + b1/D[cm])/1000 [kg/s], only segments > 1.1 cm (as PWDB)
%     terminals: 3-element Windkessel, R1 = Z0 = rho*c0/A0 at the outlet,
%                R2 = R_total - R1, C = peripheral compliance
%   Simplification: static pressure continuity at the junctions.
%
%   Numerics (staggered finite volumes, explicit, fixed step):
%     cells : area A (mass balance, exact conservation); faces: flow Q
%     symplectic (forward) Euler: Q faces first, then A cells with the new Q
%     friction semi-implicit; nodes: one linear equation each (no iterations)
%     Windkessels: forward Euler; wall viscosity: implicit diffusion on Q
%     (tridiagonal, factorised here once)
%   Stability: CFL = (|u|+c)*dt/dx <= 0.9. Segments too short for one cell
%   at this dt are merged (connectors) or removed (terminals, Windkessel
%   moved to the inlet node).
%
%   All fields of S are numeric (usable as a non-tunable parameter of a
%   MATLAB Function block).

if nargin < 5 || isempty(region), region = 'full'; end
if nargin < 6, stubF = []; end

mmHg = 133.322;
G = csvread(geoFile, 1, 0);
segIn = G(:,2); segOut = G(:,3); Ls = G(:,4);
rin = G(:,5); rout = G(:,6); Cp = G(:,7); Rp = G(:,8);
nSeg = size(G, 1);
rho = prm.rho; mu = prm.mu;
CFL_MAX = 0.9;  U_MAX = 1.5;                 % [-], [m/s] margin for blood velocity

% ---- region: roots and segments ----
F = zeros(nSeg,1);
switch lower(region)
    case 'full',      roots = 1;
    case 'upper',     roots = [4 5 15 19];
        if numel(stubF) ~= 4, error('charlton1d_setup: stubF must have 4 values for ''upper''.'); end
        F(roots) = min(max(stubF(:), 0), 1);
    case 'abdominal', roots = [26 29 34 36 38 40];
    case 'lower',     roots = [42 43];
    otherwise, error('charlton1d_setup: unknown region ''%s''.', region);
end
inReg = false(nSeg,1);  rootNodes = [];  reached = [];
for s = roots
    if 1 - F(s) < 1e-3                       % fully physical: start from its outlet
        rootNodes(end+1) = segOut(s); reached(end+1) = segOut(s); %#ok<AGROW>
    else                                     % keep only the virtual remainder
        rin(s) = rin(s) + F(s)*(rout(s) - rin(s));
        Ls(s)  = (1 - F(s))*Ls(s);
        inReg(s) = true;  rootNodes(end+1) = segIn(s); reached(end+1) = segOut(s); %#ok<AGROW>
    end
end
changed = true;                              % everything downstream of the roots
while changed
    changed = false;
    for s = 1:nSeg
        if ~inReg(s) && ~any(roots == s) && any(reached == segIn(s))
            inReg(s) = true;  reached(end+1) = segOut(s);  changed = true; %#ok<AGROW>
        end
    end
end

% ---- node merging: all root nodes -> one; short connectors -> one node ----
Kst = ones(nSeg,1);                          % stiffness factor per segment (default 1)
if isfield(prm, 'Kstiff'), Kst = prm.Kstiff(:); end
cMax  = sqrt(Kst*2/3/rho*0.1.*(prm.k1*exp(prm.k2*100*min(rin,rout)) + prm.k3));
dxMin = (cMax + U_MAX)*dt/CFL_MAX;          % smallest stable cell per segment
tooShort = inReg & Ls < dxMin;
short  = tooShort & Rp == 0;                 % connectors: merged
tShort = tooShort & Rp > 0;                  % short terminals: removed, Windkessel
                                             % moved to the inlet node
par = 1:max([segIn; segOut]);               % union-find on geometry node ids
for k = 2:numel(rootNodes), par = unite(par, rootNodes(1), rootNodes(k)); end
for s = find(short)', par = unite(par, segIn(s), segOut(s)); end
for n = 1:numel(par), par(n) = rootOf(par, n); end
segIn = par(segIn); segOut = par(segOut);
segIn = segIn(:); segOut = segOut(:);
keep = find(inReg & ~tooShort);  nK = numel(keep);

% ---- cells and faces ----
N = zeros(nSeg,1);
N(keep) = max(1, min(ceil(Ls(keep)/dxMax), floor(Ls(keep)./dxMin(keep))));
nC = sum(N);  nF = sum(N(keep) + 1);
cSeg = zeros(nC,1); A0 = cSeg; beta = cSeg; dxc = cSeg;
fSeg = zeros(nF,1); fdx = fSeg; fnu = fSeg;
cFirst = zeros(nSeg,1); fFirst = cFirst;
kc = 0; kf = 0;
for s = keep'
    n   = N(s);  dx = Ls(s)/n;
    xc  = ((1:n)' - 0.5)/n;                  % cell centres
    xf  = (0:n)'/n;                          % faces
    rc  = rin(s) + (rout(s) - rin(s))*xc;
    rf  = rin(s) + (rout(s) - rin(s))*xf;
    ic  = kc + (1:n)';  iff = kf + (1:n+1)';
    A0(ic)   = pi*rc.^2;
    Eh       = Kst(s)*0.1*(prm.k1*exp(prm.k2*100*rc) + prm.k3).*rc;   % [Pa m]
    beta(ic) = 4/3*sqrt(pi)*Eh;
    dxc(ic)  = dx;  cSeg(ic) = s;
    visc     = G(s,4) > 0.011;               % PWDB rule, on the original length
    gf       = visc*(prm.b0 + prm.b1./(200*rf))/1000;            % [kg/s]
    fdx(iff) = dx;  fSeg(iff) = s;
    fnu(iff) = gf./(rho*sqrt(pi*rf.^2));                         % [m^2/s]
    cFirst(s) = ic(1);  fFirst(s) = iff(1);
    kc = ic(end);  kf = iff(end);
end
sqA0 = sqrt(A0);
cK = sqrt(beta./(2*rho*A0));  c0 = cK.*A0.^0.25;                % wave speed at A0
cfl = max((c0 + U_MAX)*dt./dxc);
if cfl > CFL_MAX
    error('charlton1d_setup: CFL = %.2f > %.2f. Reduce dt or increase dxMax.', cfl, CFL_MAX);
end

% ---- face <-> cell connectivity: face j has left cell j-1, right cell j ----
cLast = cFirst + N - 1;  fLast = fFirst + N;
fL = zeros(nF,1); fR = fL;                   % 0 = network node side
for s = keep'
    iff = fFirst(s):fLast(s);  ic = cFirst(s):cLast(s);
    fL(iff(2:end))   = ic;
    fR(iff(1:end-1)) = ic;
end
fInt = find(fL > 0 & fR > 0);                % interior faces
cLF = zeros(nC,1); cRF = cLF;                % left/right face of each cell
cLF(fR(fR > 0)) = find(fR > 0);
cRF(fL(fL > 0)) = find(fL > 0);

% ---- end faces and network nodes ----
eF   = [fFirst(keep); fLast(keep)];          % end faces
eS   = [-ones(nK,1); ones(nK,1)];            % -1 inlet face, +1 outlet face
eC   = [cFirst(keep); cLast(keep)];          % adjacent cell
eGeo = [segIn(keep); segOut(keep)];
[geoIds, ~, eNode] = unique(eGeo);  eNode = eNode(:);  geoIds = geoIds(:);
nNodes = numel(geoIds);
root  = find(geoIds == par(rootNodes(1)));
if isempty(root), error('charlton1d_setup: root node not found.'); end
tSeg  = find(inReg & Rp > 0);
tGeo  = segOut(tSeg);  tGeo(tShort(tSeg)) = segIn(tSeg(tShort(tSeg)));
[~, tNode] = ismember(tGeo, geoIds);
if any(tNode == 0), error('charlton1d_setup: Windkessel without node.'); end
deg = accumarray(eNode, 1, [nNodes 1]);
lone = find(deg == 1); lone(ismember(lone, [root; tNode])) = [];
if ~isempty(lone), error('charlton1d_setup: dangling nodes %s.', mat2str(geoIds(lone)')); end

% ---- Windkessels: R1 = characteristic impedance at the outlet ----
rO = rout(tSeg);
cO = sqrt(Kst(tSeg)*2/3/rho*0.1.*(prm.k1*exp(prm.k2*100*rO) + prm.k3));
R1 = rho*cO./(pi*rO.^2);
R2 = Rp(tSeg) - R1;
if any(R2 <= 0), error('charlton1d_setup: R2 <= 0.'); end

% ---- implicit viscous step: (I - dt*nu*D2) Q = Q*, tridiagonal ----
% interior faces: [-w, 1+2w, -w]; end faces: identity (decoupled segments)
tdL = zeros(nF,1); tdU = tdL; tdD = ones(nF,1);
w = dt*fnu(fInt)./fdx(fInt).^2;
tdL(fInt) = -w;  tdU(fInt) = -w;  tdD(fInt) = 1 + 2*w;
% Thomas factorisation (constant matrix): store c' and 1/m
tc = zeros(nF,1); tm = zeros(nF,1);
tm(1) = 1/tdD(1);  tc(1) = tdU(1)*tm(1);
for i = 2:nF
    tm(i) = 1/(tdD(i) - tdL(i)*tc(i-1));
    tc(i) = tdU(i)*tm(i);
end

alpha = prm.alpha;
S = struct('dt', dt, 'rho', rho, 'alpha', alpha, 'Kr', 2*pi*mu*alpha/(rho*(alpha - 1)), ...
    'Pref', prm.Pref*mmHg, 'Pout', prm.Pout*mmHg, ...
    'nC', nC, 'nF', nF, 'cSeg', cSeg, 'A0', A0, 'sqA0', sqA0, 'beta', beta, ...
    'cK', cK, 'c0', c0, 'dxc', dxc, ...
    'fSeg', fSeg, 'fdx', fdx, 'fL', fL, 'fR', fR, 'fInt', fInt, 'cLF', cLF, 'cRF', cRF, ...
    'cFirst', cFirst, 'cLast', cLast, 'fFirst', fFirst, 'fLast', fLast, ...
    'eF', eF, 'eS', eS, 'eC', eC, 'eNode', eNode, 'nNodes', nNodes, ...
    'geoIds', geoIds, 'root', root, 'tSeg', tSeg, 'tNode', tNode, ...
    'R1', R1, 'R2', R2, 'CT', Cp(tSeg), 'tdL', tdL, 'tc', tc, 'tm', tm, ...
    'N', N, 'merged', find(short), 'removedTerminal', find(tShort), 'cfl', cfl);

fprintf(['1D network (%s): %d segments (merged %s, terminal removed %s), %d cells, ' ...
         '%d nodes, %d Windkessels, dt = %g s, dx <= %g cm, CFL <= %.2f\n'], ...
        region, nK, mat2str(find(short)'), mat2str(find(tShort)'), nC, nNodes, ...
        numel(tSeg), dt, 100*dxMax, cfl);
end

function a = rootOf(par, a)
while par(a) ~= a, a = par(a); end
end

function par = unite(par, a, b)
a = rootOf(par, a);  b = rootOf(par, b);
par(max(a,b)) = min(a,b);
end
