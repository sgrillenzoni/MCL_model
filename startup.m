%% startup.m
% Adds the repository folders to the MATLAB path.
% Run once from the repository root (or open MATLAB in this folder, where
% it runs automatically) before opening any model.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'params'), fullfile(root, 'solver1d'), fullfile(root, 'blocks'), ...
        fullfile(root, 'models'), fullfile(root, 'models', 'tests'), ...
        fullfile(root, 'calibration'), fullfile(root, 'validation'), ...
        fullfile(root, 'data'), fullfile(root, 'external', 'pwdb'));
fprintf('MCL_model: paths added (%s)\n', root);
clear root
