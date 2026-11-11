%% compare_pwdb_1d.m
% Compare the 1D solver (run_charlton_1d.m) with the PWDB reference
% simulation of the same subject (Charlton et al. 2019, baseline 25 y),
% at the same sites, over the last beat.
%
% Reference data: Zenodo record 3472909 ("PWDB: baseline subjects aged
% 25 to 75"), simulation sim_1 (25 y, HR 72.9 bpm), Nektar1D history
% files sim_1_<segment>.his. Columns: t [s], P [Pa], Pe [Pa], U [m/s],
% Q [m^3/s], A [m^2], history point. The files needed are downloaded once
% into the folder PWDB_DIR (a few MB each).
%
% Run AFTER run_charlton_1d (needs res, sites, CH_T, CH_LEN in the workspace).

PWDB_DIR = fullfile(fileparts(which('pwdb_geo0001.csv')), 'pwdb_sim1');   % data/pwdb_sim1 (not tracked, see README)
URL      = 'https://zenodo.org/record/3472909/files/';
SITES    = 1:6;                      % AorticRoot, Carotid, Brachial, Radial, Femoral, AbdAorta
mmHg     = 133.322;

if ~exist(PWDB_DIR, 'dir'), mkdir(PWDB_DIR); end
T  = CH_T;
ph = (0:1e-3:T)';                    % common phase grid [s]

% our last complete beat, as phase within the cycle
nb  = floor(res.t(end)/T);
iT  = res.t > (nb-1)*T & res.t <= nb*T;
phO = mod(res.t(iT), T);  [phO, o] = sort(phO);

nS = numel(SITES);
ref = struct('P', zeros(numel(ph), nS), 'Q', zeros(numel(ph), nS));
our = ref;
fprintf('\n%-11s | %-21s | %-21s | %-24s | %s\n', '', 'SBP / DBP / MAP [mmHg]', ...
        'PWDB SBP / DBP / MAP', 'Qmean  back%  (1D | PWDB)', 'RMSE P [mmHg]  Q [mL/s]');
for k = 1:nS
    j = SITES(k);  seg = sites{j,2};  pos = sites{j,3};
    f = fullfile(PWDB_DIR, sprintf('sim_1_%d.his', seg));
    if ~exist(f, 'file')
        fprintf('Downloading %s ...\n', f);
        if exist('websave', 'file'), websave(f, [URL sprintf('sim_1_%d.his', seg)]);
        else, urlwrite([URL sprintf('sim_1_%d.his', seg)], f); end
    end
    % history points and their positions
    txt = fileread(f);
    xs  = regexp(txt, '# Point \d+: x = ([\d.]+)', 'tokens');
    xs  = cellfun(@(c) str2double(c{1}), xs);
    [~, pt] = min(abs(xs - pos*CH_LEN(seg)/100));
    fid = fopen(f);  D = textscan(fid, '%f %f %f %f %f %f %f', 'CommentStyle', '#');  fclose(fid);
    D = [D{:}];  D = D(D(:,7) == pt, :);
    % PWDB last beat, as phase within the cycle
    tR = D(:,1);  iR = tR > tR(end) - T;
    phR = mod(tR(iR), T);  [phR, oR] = sort(phR);
    PR = D(iR, 2)/mmHg;  QR = 1e6*D(iR, 5);
    [phR, u] = unique(phR);  PR = PR(oR(u));  QR = QR(oR(u));
    ref.P(:,k) = interp1(phR, PR, ph, 'linear', 'extrap');
    ref.Q(:,k) = interp1(phR, QR, ph, 'linear', 'extrap');
    PO = res.P(iT, j);  QO = res.Q(iT, j);
    [phU, u] = unique(phO);  PO = PO(o(u));  QO = QO(o(u));
    our.P(:,k) = interp1(phU, PO, ph, 'linear', 'extrap');
    our.Q(:,k) = interp1(phU, QO, ph, 'linear', 'extrap');

    m = @(y) trapz(ph, y)/T;  bf = @(y) 100*max(0, -min(y))/max(y);
    fprintf('%-11s | %5.1f / %5.1f / %5.1f   | %5.1f / %5.1f / %5.1f   | %6.2f %4.0f | %6.2f %4.0f   | %5.1f   %6.2f\n', ...
        sites{j,1}, max(our.P(:,k)), min(our.P(:,k)), m(our.P(:,k)), ...
        max(ref.P(:,k)), min(ref.P(:,k)), m(ref.P(:,k)), ...
        m(our.Q(:,k)), bf(our.Q(:,k)), m(ref.Q(:,k)), bf(ref.Q(:,k)), ...
        sqrt(mean((our.P(:,k) - ref.P(:,k)).^2)), sqrt(mean((our.Q(:,k) - ref.Q(:,k)).^2)));
end

figure('Name', '1D solver vs PWDB (25 y), last beat');
for k = 1:nS
    subplot(2, nS, k);
    plot(ph, ref.P(:,k), 'k', ph, our.P(:,k), 'r'); grid on
    title(sites{SITES(k),1}, 'Interpreter', 'none');
    if k == 1, ylabel('P [mmHg]'); legend('PWDB', '1D solver'); end
    subplot(2, nS, nS + k);
    plot(ph, ref.Q(:,k), 'k', ph, our.Q(:,k), 'r', ph, 0*ph, 'k:'); grid on
    xlabel('t [s]');  if k == 1, ylabel('Q [mL/s]'); end
end
