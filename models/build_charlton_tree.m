function sub = build_charlton_tree(parent, name, region, maxCell)
%BUILD_CHARLTON_TREE  Build the virtual arterial tree of one Windkessel
%   chamber from the PWDB network (Charlton et al. 2019, 116 segments) as a
%   Simscape subsystem with ONE physical port 'in' (the chamber node).
%
%   sub = build_charlton_tree(parent, name, region)
%   sub = build_charlton_tree(parent, name, region, maxCell)
%
%   parent  : system where the subsystem is created (e.g. gcs)
%   name    : name of the new subsystem (e.g. 'UpperBody_tree')
%   region  : 'upper' | 'abdominal' | 'lower'
%   maxCell : maximum cell length [cm] (default 3). Each segment is split
%             into N = ceil(length/maxCell) identical cells.
%
%   Unlike build_vascular_tree.m, the network is built NODE BY NODE from the
%   inlet/outlet node numbers of the PWDB geometry, so loops are allowed
%   (circle of Willis, palmar arches).
%
%   Each cell:   --[R/N]--[L/N]--o--          o = cell node
%                                |
%                             [Rv]   wall viscosity, Rv = KTAU_V*TAU / C_cell
%                                |
%                             [C/N]
%                                |
%                               ref
%   Each terminal segment ends with the PWDB 3-element Windkessel:
%       node --[R1]--+--[R2]-- P_venous
%                    |
%                   [CT]
%                    |
%                   ref
%
%   All values are EXPRESSIONS of workspace variables (CH_R, CH_L, CH_C,
%   CH_TAU, CH_R1, CH_R2, CH_CT, P_venous, KTAU_V, KC_VIRTUAL,
%   F_STUB_CAROTID_L) -> run parameters_charlton.m before building and
%   before simulating (model callbacks).
%
%   Block names: R_s<seg>_<cell>, L_s<seg>_<cell>, ... e.g. R_s5_1 is the
%   first cell of segment 5 (R common carotid).
%
%   Regions (what is physical in the phantom, or aorta, is NOT included):
%     upper     : the virtual remainder (1 - F_STUB_...) of the four phantom
%                 outlets: 4 (R subclavian A), 5 (R common carotid),
%                 15 (L common carotid), 19 (L subclavian A, before the
%                 L vertebral) + everything downstream (arms, hands, head).
%                 The innominate (3) is fully physical and not included.
%     abdominal : roots 26 (intercostals), 29 (celiac), 34 (sup. mesenteric),
%                 36, 38 (renals), 40 (inf. mesenteric) + downstream
%     lower     : roots 42, 43 (common iliacs) + downstream
%   The inlet nodes of all roots are joined into the single port 'in'.
%
%   Uses only Simscape Foundation electrical blocks available in R2018.
%   Example:
%     parameters_charlton;
%     build_charlton_tree(gcs, 'UpperBody_tree', 'upper');

if nargin < 4 || isempty(maxCell), maxCell = 3; end

% ---------------- network data from the workspace ----------------
need = {'CH_IN', 'CH_OUT', 'CH_LEN', 'CH_TERM'};
for k = 1:numel(need)
    if ~evalin('base', sprintf('exist(''%s'', ''var'')', need{k}))
        error('build_charlton_tree: %s not found. Run parameters_charlton.m first.', need{k});
    end
end
IN   = evalin('base', 'CH_IN');
OUT  = evalin('base', 'CH_OUT');
LEN  = evalin('base', 'CH_LEN');
TERM = evalin('base', 'CH_TERM');

% ---------------- region ----------------
% Roots of the upper body = the four phantom outlets. Only the virtual
% remainder (1 - F_STUB_...) of each root segment is in the tree; a root
% that is fully physical (F = 1) is skipped and its outlet node becomes
% the port 'in' (e.g. L subclavian A fully physical -> vertebral and
% subclavian B start at 'in').
frac = containers.Map('KeyType', 'double', 'ValueType', 'any');
switch lower(region)
    case 'upper'
        roots = [4 5 15 19];
        frac(4)  = 'F_STUB_SUBCLAVIAN_R';  frac(5)  = 'F_STUB_CAROTID_R';
        frac(15) = 'F_STUB_CAROTID_L';     frac(19) = 'F_STUB_SUBCLAVIAN_L';
    case 'abdominal', roots = [26 29 34 36 38 40];
    case 'lower',     roots = [42 43];
    otherwise
        error('build_charlton_tree: region must be ''upper'', ''abdominal'' or ''lower''.');
end
fracNum = containers.Map('KeyType', 'double', 'ValueType', 'double');   % virtual fraction
inNodes = [];  keep = [];
for s = roots
    f = 0;
    if isKey(frac, s)
        if ~evalin('base', sprintf('exist(''%s'', ''var'')', frac(s)))
            error('build_charlton_tree: %s not found. Run parameters_charlton.m first.', frac(s));
        end
        f = evalin('base', frac(s));
        frac(s) = sprintf('(1-%s)', frac(s));
    end
    if 1 - f < 1e-3                          % fully physical: skip, start from its outlet
        inNodes(end+1) = OUT(s); %#ok<AGROW>
    else
        fracNum(s) = 1 - f;  keep(end+1) = s;  inNodes(end+1) = IN(s); %#ok<AGROW>
    end
end
roots = keep;

% segments of the region: kept roots + everything reachable downstream
segs = roots;  reached = [OUT(roots) inNodes(~ismember(inNodes, IN(roots)))];
for s = 1:numel(IN)                          % children of skipped roots
    if ~any(segs == s) && any(inNodes == IN(s)) && ~any(IN(roots) == IN(s))
        segs(end+1) = s; reached(end+1) = OUT(s); %#ok<AGROW>
    end
end
changed = true;
while changed
    changed = false;
    for s = 1:numel(IN)
        if ~any(segs == s) && any(reached == IN(s))
            segs(end+1) = s;  reached(end+1) = OUT(s);  changed = true; %#ok<AGROW>
        end
    end
end

% ---------------- libraries ----------------
load_system('fl_lib'); load_system('nesl_utility');
LIB.R    = 'fl_lib/Electrical/Electrical Elements/Resistor';
LIB.L    = 'fl_lib/Electrical/Electrical Elements/Inductor';
LIB.C    = 'fl_lib/Electrical/Electrical Elements/Capacitor';
LIB.Ref  = 'fl_lib/Electrical/Electrical Elements/Electrical Reference';
LIB.Vdc  = 'fl_lib/Electrical/Electrical Sources/DC Voltage Source';
LIB.Port = 'nesl_utility/Connection Port';

% ---------------- create the subsystem ----------------
root = strtok(parent, '/'); load_system(root);
sub = [parent '/' name];
if getSimulinkBlockHandle(sub) ~= -1
    error('build_charlton_tree: "%s" already exists. Delete it or choose another name.', sub);
end
add_block('built-in/Subsystem', sub, 'Position', [100 100 260 160]);
set_param(sub, 'Description', sprintf(['Virtual arterial tree (%s), PWDB baseline 25 y ' ...
    '(Charlton 2019), cells <= %g cm. Generated by build_charlton_tree.m'], region, maxCell));
pIn = add_block(LIB.Port, [sub '/in'], 'Position', [20 95 50 115], 'Side', 'left');

% node map: node number -> one physical port handle on that node
node = containers.Map('KeyType', 'double', 'ValueType', 'double');
hIn  = portOf(pIn);
for n0 = inNodes, node(n0) = hIn; end

% ---------------- build, node by node ----------------
CELL_W = 300; ROW_H = 260; row = 0; nCells = 0; nBlk = 0;
built = false(1, numel(IN));
queue = unique(inNodes, 'stable');
while ~isempty(queue)
    n = queue(1); queue(1) = [];
    for s = segs
        if built(s) || IN(s) ~= n, continue; end
        built(s) = true;
        y = 100 + row*ROW_H; row = row + 1;
        addSegment(s, y);
        if ~any(queue == OUT(s)), queue(end+1) = OUT(s); end %#ok<AGROW>
    end
end
if ~all(built(segs))
    error('build_charlton_tree: segments %s not reached.', mat2str(segs(~built(segs))));
end
fprintf('Created %s: %d segments, %d cells, %d terminals, %d blocks (cells <= %g cm).\n', ...
        sub, numel(segs), nCells, nnz(TERM(segs)), nBlk, maxCell);

% =====================================================================
    function addSegment(s, y)
        fr = '1'; if isKey(frac, s), fr = frac(s); end
        fn = 1; if isKey(fracNum, s), fn = fracNum(s); end
        N  = max(1, ceil(fn*LEN(s)/maxCell));
        cExpr = sprintf('(KC_VIRTUAL*%s*CH_C(%d)/%d)', fr, s, N);
        h = node(IN(s)); x = 80;
        for k = 1:N
            tag = sprintf('s%d_%d', s, k);
            bR  = blk(LIB.R,   ['R_' tag],  [x 0 x+40 30]       + [0 y 0 y], 'R', sprintf('%s*CH_R(%d)/%d', fr, s, N));
            bL  = blk(LIB.L,   ['L_' tag],  [x+70 0 x+110 30]   + [0 y 0 y], 'l', sprintf('%s*CH_L(%d)/%d', fr, s, N));
            bRv = blk(LIB.R,   ['Rv_' tag], [x+150 60 x+190 90] + [0 y 0 y], 'R', ...
                      sprintf('max(KTAU_V*CH_TAU(%d),1e-9)/%s', s, cExpr));
            bC  = blk(LIB.C,   ['C_' tag],  [x+150 120 x+190 150] + [0 y 0 y], 'c', cExpr);
            bG  = blk(LIB.Ref, ['Ref_' tag],[x+150 180 x+190 210] + [0 y 0 y]);
            pR = ph(bR); pL = ph(bL); pRv = ph(bRv); pC = ph(bC); pG = ph(bG);
            ln(h, pR.LConn(1)); ln(pR.RConn(1), pL.LConn(1));
            ln(pL.RConn(1), pRv.LConn(1)); ln(pRv.RConn(1), pC.LConn(1)); ln(pC.RConn(1), pG.LConn(1));
            h = pL.RConn(1); x = x + CELL_W; nCells = nCells + 1;
        end
        % outlet node: join an existing node (loops) or create it
        if isKey(node, OUT(s))
            ln(h, node(OUT(s)));
        else
            node(OUT(s)) = h;
        end
        if TERM(s)                                   % PWDB 3-element Windkessel
            t   = sprintf('s%d', s);
            bR1 = blk(LIB.R,   ['R1_' t],   [x 0 x+40 30]        + [0 y 0 y], 'R', sprintf('CH_R1(%d)', s));
            bR2 = blk(LIB.R,   ['R2_' t],   [x+80 0 x+120 30]    + [0 y 0 y], 'R', sprintf('CH_R2(%d)', s));
            bV  = blk(LIB.Vdc, ['Pv_' t],   [x+160 0 x+200 30]   + [0 y 0 y], 'v0', 'P_venous');
            bGv = blk(LIB.Ref, ['RefV_' t], [x+230 0 x+270 30]   + [0 y 0 y]);
            bCT = blk(LIB.C,   ['CT_' t],   [x+60 120 x+100 150] + [0 y 0 y], 'c', sprintf('KC_VIRTUAL*CH_CT(%d)', s));
            bGc = blk(LIB.Ref, ['RefCT_' t],[x+60 180 x+100 210] + [0 y 0 y]);
            p1 = ph(bR1); p2 = ph(bR2); pV = ph(bV); pGv = ph(bGv); pCT = ph(bCT); pGc = ph(bGc);
            ln(node(OUT(s)), p1.LConn(1));
            ln(p1.RConn(1), p2.LConn(1)); ln(p2.RConn(1), pV.LConn(1)); ln(pV.RConn(1), pGv.LConn(1));
            ln(p1.RConn(1), pCT.LConn(1)); ln(pCT.RConn(1), pGc.LConn(1));
        end
    end

    function b = blk(lib, nm, pos, varargin)
        b = add_block(lib, [sub '/' nm], 'Position', pos, varargin{:});
        nBlk = nBlk + 1;
    end

    function ln(h1, h2)
        add_line(sub, h1, h2, 'autorouting', 'on');
    end
end

% =====================================================================
function p = ph(b)
p = get_param(b, 'PortHandles');
end

function h = portOf(blk)
p = get_param(blk, 'PortHandles');
h = [p.LConn p.RConn]; h = h(1);
end
