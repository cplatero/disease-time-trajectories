% setup.m  Adds src/ and experiments/ to the MATLAB path.
% Run it once per session from the root of the repository, before calling
% any experiment or the tutorial directly (run_all does it by itself).
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'src'), fullfile(root, 'experiments'), fullfile(root, 'tutorial'));
clear root
