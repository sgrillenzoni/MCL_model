function sub = build_segment_cells(parent, name, Rvar, Lvar, Cvar, N, pos)
%BUILD_SEGMENT_CELLS  Build a Simscape subsystem that represents one arterial
%   segment as N identical R-L-C cells in series (lumped transmission line).
%
%   sub = build_segment_cells(parent, name, Rvar, Lvar, Cvar, N)
%   sub = build_segment_cells(parent, name, Rvar, Lvar, Cvar, N, pos)
%
%   parent : system where the subsystem is created, e.g. 'mockloop_numerical/UpperBody'
%   name   : name of the new subsystem, e.g. 'subclavianB_right_10cells'
%   Rvar, Lvar, Cvar : NAMES of the workspace variables of the whole segment
%            (as in parameters_stergiopulos.m), e.g. 'R_subclavianB_right'
%   N      : number of cells
%   pos    : (optional) [left top right bottom] of the subsystem in parent
%
%   Why: a long vessel modelled as ONE R-L-C cell behaves like a single
%   mass-spring oscillator with a resonance at f = 1/(2*pi*sqrt(L*C)). For the
%   subclavian B (42 cm) and the femoral (44 cm) this is 2.2-2.5 Hz, on the
%   second harmonic of the heartbeat at 75 bpm, and it produces large
%   artificial flow oscillations. Splitting the segment into N cells keeps
%   the totals (sum of R, L and C unchanged, so mean flows and pressures do
%   not change) while each cell resonates at about N times higher frequency,
%   outside the band of the heartbeat. Rule of thumb: cells of 4-5 cm at most.
%
%   Each cell, as in the existing network:  --[R/N]--[L/N]--o--
%                                                           |
%                                                         [C/N]
%                                                           |
%                                                          ref
%   The subsystem has two physical ports: 'in' (left, proximal end) and
%   'out' (right, distal end). Values are written as expressions of the
%   workspace variables (e.g. 'R_subclavianB_right/10'), so they stay linked
%   to the parameter script.
%
%   Uses only Simscape Foundation electrical blocks available in R2018.
%
%   Example:
%     build_segment_cells('mockloop_numerical/UpperBody', 'subclavianB_right_10cells', ...
%         'R_subclavianB_right', 'L_subclavianB_right', 'C_subclavianB_right', 10);

if nargin < 7 || isempty(pos), pos = [100 100 220 160]; end

% Libraries
load_system('fl_lib');          % Simscape Foundation Library
load_system('nesl_utility');    % Simscape utilities (Connection Port)
libR   = 'fl_lib/Electrical/Electrical Elements/Resistor';
libL   = 'fl_lib/Electrical/Electrical Elements/Inductor';
libC   = 'fl_lib/Electrical/Electrical Elements/Capacitor';
libRef = 'fl_lib/Electrical/Electrical Elements/Electrical Reference';
libPort= 'nesl_utility/Connection Port';

% Root model must be loaded
root = strtok(parent, '/');
load_system(root);

% Create an empty subsystem
sub = [parent '/' name];
if getSimulinkBlockHandle(sub) ~= -1
    error('build_segment_cells: block "%s" already exists. Delete it or choose another name.', sub);
end
add_block('built-in/Subsystem', sub, 'Position', pos);
set_param(sub, 'Description', sprintf( ...
    'Segment %s / %s / %s split into %d R-L-C cells (each = total/%d).', ...
    Rvar, Lvar, Cvar, N, N));

% Ports
dx = 220;                                   % horizontal spacing of the cells
pIn  = add_block(libPort, [sub '/in'],  'Position', [20 95 50 115], 'Side', 'left');
pOut = add_block(libPort, [sub '/out'], 'Position', [80 + N*dx 95 110 + N*dx 115], 'Side', 'right');

prev = portOf(pIn);                         % node to connect the next cell to
for k = 1:N
    x0 = 80 + (k-1)*dx;
    hR = add_block(libR, sprintf('%s/R%d', sub, k), 'Position', [x0      90 x0+40  120], ...
                   'R', sprintf('%s/%d', Rvar, N));
    hL = add_block(libL, sprintf('%s/L%d', sub, k), 'Position', [x0+70   90 x0+110 120], ...
                   'l', sprintf('%s/%d', Lvar, N));
    hC = add_block(libC, sprintf('%s/C%d', sub, k), 'Position', [x0+150 150 x0+180 190], ...
                   'Orientation', 'down', 'c', sprintf('%s/%d', Cvar, N));
    hG = add_block(libRef, sprintf('%s/Ref%d', sub, k), 'Position', [x0+150 230 x0+180 260]);

    pR = get_param(hR, 'PortHandles');
    pL = get_param(hL, 'PortHandles');
    pC = get_param(hC, 'PortHandles');
    pG = get_param(hG, 'PortHandles');

    add_line(sub, prev,       pR.LConn(1), 'autorouting', 'on');   % previous node -> R+
    add_line(sub, pR.RConn(1), pL.LConn(1), 'autorouting', 'on');  % R- -> L+
    add_line(sub, pL.RConn(1), pC.LConn(1), 'autorouting', 'on');  % L- -> C+ (cell node)
    add_line(sub, pC.RConn(1), pG.LConn(1), 'autorouting', 'on');  % C- -> reference

    prev = pL.RConn(1);                     % the cell node feeds the next cell
end
add_line(sub, prev, portOf(pOut), 'autorouting', 'on');

fprintf('Created %s: %d cells of %s/%d, %s/%d, %s/%d.\n', sub, N, Rvar, N, Lvar, N, Cvar, N);
end

function h = portOf(blk)
% The single physical port of a Connection Port block.
ph = get_param(blk, 'PortHandles');
h  = [ph.LConn ph.RConn];
h  = h(1);
end
