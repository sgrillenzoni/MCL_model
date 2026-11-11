function sub = build_vascular_tree(parent, name, region, maxCell)
%BUILD_VASCULAR_TREE  Build automatically the virtual arterial tree of one
%   Windkessel chamber (Stergiopulos et al. 1992 network, phantom topology) as a
%   Simscape subsystem with ONE physical port 'in' (the chamber node).
%
%   sub = build_vascular_tree(parent, name, region)
%   sub = build_vascular_tree(parent, name, region, maxCell)
%
%   parent  : system where the subsystem is created (e.g. gcs)
%   name    : name of the new subsystem (e.g. 'UpperBody_tree')
%   region  : 'upper' | 'abdominal' | 'lower'
%   maxCell : maximum cell length [cm] (default 3). Every vessel is split into
%             N = ceil(length/maxCell) identical R-L-C cells, so that no cell
%             resonates inside the band of the heartbeat.
%
%   Each cell:   --[R/N]--[L/N]--o--          o = cell node
%                                |
%                             [Rv]   wall viscosity (Voigt), Rv = TAU_V / C_cell
%                                |
%                             [C/N]
%                                |
%                               ref
%   Each terminal vessel ends with its peripheral bed:
%       node --[RT]-- P_venous      and      node --[CT]-- ref
%
%   All values are EXPRESSIONS of workspace variables, so they can be changed
%   without rebuilding the block:
%     R_<vessel>, L_<vessel>, C_<vessel>, RT_<vessel>, CT_<vessel>, P_venous
%                      -> parameters_stergiopulos.m
%     TAU_V, KC_VIRTUAL, KRT_HEAD, KCT_LEGS, F_STUB_CAROTID_L
%                      -> parameters_vascular_tuning.m
%   Run both scripts before building and before simulating (model callbacks).
%
%   Use: place the generated subsystem INSIDE the existing UpperBody /
%   Abdominal / LowerBody block, in place of the old network, and connect its
%   'in' port to the node where the old network started (the node fed by the
%   flow from the twin in MODE 1 or by the controlled flow source in MODE 2).
%   The MODE switching of the existing block is not touched.
%
%   Topology (what is physical in the phantom is NOT included):
%     upper     : R subclavian A, R common carotid, rest of L common carotid,
%                 L vertebral, L subclavian B, with their distal trees
%                 (innominate and L subclavian A are in the phantom)
%     abdominal : intercostals, celiac tree, superior mesenteric, renals,
%                 inferior mesenteric (all from the renal outlets of the phantom)
%     lower     : common iliacs with their distal trees
%
%   Uses only Simscape Foundation electrical blocks available in R2018.
%   Example:
%     parameters_stergiopulos; parameters_vascular_tuning;
%     build_vascular_tree(gcs, 'UpperBody_tree', 'upper');

if nargin < 4 || isempty(maxCell), maxCell = 3; end

% ---------------- libraries ----------------
load_system('fl_lib'); load_system('nesl_utility');
LIB.R    = 'fl_lib/Electrical/Electrical Elements/Resistor';
LIB.L    = 'fl_lib/Electrical/Electrical Elements/Inductor';
LIB.C    = 'fl_lib/Electrical/Electrical Elements/Capacitor';
LIB.Ref  = 'fl_lib/Electrical/Electrical Elements/Electrical Reference';
LIB.Vdc  = 'fl_lib/Electrical/Electrical Sources/DC Voltage Source';
LIB.Port = 'nesl_utility/Connection Port';

% ---------------- vessel lengths [cm], Stergiopulos Table 1 ----------------
len = struct( ...
 'subclavianA_right',3.4, 'commonCarotid_right',17.7, 'vertebral_right',14.8, ...
 'subclavianB_right',42.2, 'radial_right',23.5, 'ulnarA_right',6.7, 'interosseous_right',7.9, ...
 'ulnarB_right',17.1, 'internalCarotid_right',17.7, 'externalCarotid_right',17.7, ...
 'commonCarotid_left',20.8, 'internalCarotid_left',17.7, 'externalCarotid_left',17.7, ...
 'vertebral_left',14.8, 'subclavianB_left',42.2, 'radial_left',23.5, 'ulnarA_left',6.7, ...
 'interosseous_left',7.9, 'ulnarB_left',17.1, 'intercostals',8, 'celiacA',1, 'celiacB',1, ...
 'hepatic',6.6, 'gastric',7.1, 'splenic',6.3, 'superiorMesenteric',5.9, 'renal_left',3.2, ...
 'renal_right',3.2, 'inferiorMesenteric',5, 'commonIliac_right',5.8, 'commonIliac_left',5.8, ...
 'externalIliac_left',14.4, 'internalIliac_left',5, 'femoral_left',44.3, 'deepFemoral_left',12.6, ...
 'posteriorTibial_left',32.1, 'anteriorTibial_left',34.3, 'externalIliac_right',14.4, ...
 'internalIliac_right',5, 'femoral_right',44.3, 'deepFemoral_right',12.6, ...
 'posteriorTibial_right',32.1, 'anteriorTibial_right',34.3);

% terminal groups for the tuning factors
HEAD = {'vertebral_right','internalCarotid_right','externalCarotid_right', ...
        'internalCarotid_left','externalCarotid_left','vertebral_left'};
LEGS = {'internalIliac_left','deepFemoral_left','posteriorTibial_left','anteriorTibial_left', ...
        'internalIliac_right','deepFemoral_right','posteriorTibial_right','anteriorTibial_right'};

% ---------------- tree of the region ----------------
roots = treeSpec(region);

% ---------------- create the subsystem ----------------
root = strtok(parent, '/'); load_system(root);
sub = [parent '/' name];
if getSimulinkBlockHandle(sub) ~= -1
    error('build_vascular_tree: "%s" already exists. Delete it or choose another name.', sub);
end
add_block('built-in/Subsystem', sub, 'Position', [100 100 260 160]);
set_param(sub, 'Description', sprintf(['Virtual arterial tree (%s), Stergiopulos 1992, ' ...
    'cells <= %g cm, wall viscosity TAU_V. Generated by build_vascular_tree.m'], region, maxCell));
pIn = add_block(LIB.Port, [sub '/in'], 'Position', [20 95 50 115], 'Side', 'left');
hIn = portOf(pIn);

CELL_W = 300; ROW_H = 260;
nextRow = 100; nBlk = 0; nCells = 0;
for r = 1:numel(roots)
    y = nextRow; nextRow = nextRow + ROW_H;
    addBranch(roots{r}, hIn, 80, y);
end
fprintf('Created %s: %d cells (%d blocks), max cell length %g cm.\n', sub, nCells, nBlk, maxCell);

% =====================================================================
    function addBranch(br, hPar, x0, y)
        nm = br{1}; kids = br{3};
        fr = br{2};                                  % fraction expression ('1' or '(1-F_...)')
        frNum = 1; if ~strcmp(fr, '1'), frNum = 0.85; end   % only for the cell count
        N = max(1, ceil(frNum * len.(nm) / maxCell));
        cExpr = sprintf('(KC_VIRTUAL*%s*C_%s/%d)', fr, nm, N);
        h = hPar; x = x0;
        for k = 1:N
            tag = sprintf('%s_%d', nm, k);
            bR  = blk(LIB.R,   ['R_' tag],  [x 0 x+40 30] + [0 y 0 y],      'R', sprintf('%s*R_%s/%d', fr, nm, N));
            bL  = blk(LIB.L,   ['L_' tag],  [x+70 0 x+110 30] + [0 y 0 y],  'l', sprintf('%s*L_%s/%d', fr, nm, N));
            bRv = blk(LIB.R,   ['Rv_' tag], [x+150 60 x+190 90] + [0 y 0 y], 'R', sprintf('max(TAU_V,1e-9)/%s', cExpr));
            bC  = blk(LIB.C,   ['C_' tag],  [x+150 120 x+190 150] + [0 y 0 y], 'c', cExpr);
            bG  = blk(LIB.Ref, ['Ref_' tag],[x+150 180 x+190 210] + [0 y 0 y]);
            pR = ph(bR); pL = ph(bL); pRv = ph(bRv); pC = ph(bC); pG = ph(bG);
            ln(h, pR.LConn(1)); ln(pR.RConn(1), pL.LConn(1));
            ln(pL.RConn(1), pRv.LConn(1)); ln(pRv.RConn(1), pC.LConn(1)); ln(pC.RConn(1), pG.LConn(1));
            h = pL.RConn(1); x = x + CELL_W; nCells = nCells + 1;
        end
        if isempty(kids)                            % terminal: peripheral bed
            [kR, kC] = factors(nm);
            bRT = blk(LIB.R,   ['RT_' nm],  [x 0 x+40 30] + [0 y 0 y],       'R', sprintf('%s*RT_%s', kR, nm));
            bV  = blk(LIB.Vdc, ['Pv_' nm],  [x+80 0 x+120 30] + [0 y 0 y],   'v0', 'P_venous');
            bGv = blk(LIB.Ref, ['RefV_' nm],[x+150 0 x+190 30] + [0 y 0 y]);
            bCT = blk(LIB.C,   ['CT_' nm],  [x 120 x+40 150] + [0 y 0 y],    'c', ...
                      sprintf('KC_VIRTUAL*%s*CT_%s/%s', kC, nm, kR));
            bGc = blk(LIB.Ref, ['RefCT_' nm],[x 180 x+40 210] + [0 y 0 y]);
            pRT = ph(bRT); pV = ph(bV); pGv = ph(bGv); pCT = ph(bCT); pGc = ph(bGc);
            ln(h, pRT.LConn(1)); ln(pRT.RConn(1), pV.LConn(1)); ln(pV.RConn(1), pGv.LConn(1));
            ln(h, pCT.LConn(1)); ln(pCT.RConn(1), pGc.LConn(1));
        else
            for c = 1:numel(kids)
                if c == 1, yc = y; else, yc = nextRow; nextRow = nextRow + ROW_H; end
                addBranch(kids{c}, h, x, yc);
            end
        end
    end

    function b = blk(lib, nm, pos, varargin)
        b = add_block(lib, [sub '/' nm], 'Position', pos, varargin{:});
        nBlk = nBlk + 1;
    end

    function ln(h1, h2)
        add_line(sub, h1, h2, 'autorouting', 'on');
    end

    function [kR, kC] = factors(nm)
        kR = '1'; kC = '1';
        if any(strcmp(nm, HEAD)), kR = 'KRT_HEAD'; end
        if any(strcmp(nm, LEGS)), kC = 'KCT_LEGS'; end
    end
end

% =====================================================================
function roots = treeSpec(region)
% Branch = {vessel name, fraction expression, {children}}
B = @(n, k) {n, '1', k};
switch lower(region)
    case 'upper'
        armR = B('subclavianA_right', {B('vertebral_right', {}), ...
               B('subclavianB_right', {B('radial_right', {}), ...
               B('ulnarA_right', {B('interosseous_right', {}), B('ulnarB_right', {})})})});
        carR = B('commonCarotid_right', {B('internalCarotid_right', {}), B('externalCarotid_right', {})});
        carL = {'commonCarotid_left', '(1-F_STUB_CAROTID_L)', ...
                {B('internalCarotid_left', {}), B('externalCarotid_left', {})}};
        verL = B('vertebral_left', {});
        armL = B('subclavianB_left', {B('radial_left', {}), ...
               B('ulnarA_left', {B('interosseous_left', {}), B('ulnarB_left', {})})});
        roots = {armR, carR, carL, verL, armL};
    case 'abdominal'
        roots = {B('intercostals', {}), ...
                 B('celiacA', {B('celiacB', {B('gastric', {}), B('splenic', {})}), B('hepatic', {})}), ...
                 B('superiorMesenteric', {}), B('renal_left', {}), B('renal_right', {}), ...
                 B('inferiorMesenteric', {})};
    case 'lower'
        legR = B('commonIliac_right', {B('externalIliac_right', {B('femoral_right', ...
               {B('posteriorTibial_right', {}), B('anteriorTibial_right', {})}), ...
               B('deepFemoral_right', {})}), B('internalIliac_right', {})});
        legL = B('commonIliac_left', {B('externalIliac_left', {B('femoral_left', ...
               {B('posteriorTibial_left', {}), B('anteriorTibial_left', {})}), ...
               B('deepFemoral_left', {})}), B('internalIliac_left', {})});
        roots = {legR, legL};
    otherwise
        error('build_vascular_tree: region must be ''upper'', ''abdominal'' or ''lower''.');
end
end

function p = ph(b)
p = get_param(b, 'PortHandles');
end

function h = portOf(blk)
p = get_param(blk, 'PortHandles');
h = [p.LConn p.RConn]; h = h(1);
end
