%RUN_MODEL_SELECTION  Ostaje li 96 % za Halpin-Tsai ako se prosiri prior za E2f?
%
%   matlab -batch "run('C:/Users/imari/Documents/StipeRad/tools/run_model_selection.m')"

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
cd(root); addpath('models'); addpath('mcmc'); addpath('tools');

d = importdata('E2podaci.dat');
x = d(:, 1)'; y = d(:, 2)';

prop = [0.5 0.05 0.05];
N = 100000;
jumpPr = 0.9;

priors = [ ...
    73  5 ; ...
    65 15 ; ...
    60 15 ; ...
    55 15 ; ...
    80  5 ];

fprintf('=== Model selection vs. prior za E2f (Em~N(3,1), sigma~N(1.5,1), jumpPr=0.9, N=%d) ===\n', N);
fprintf('%9s | %-8s %-8s %-8s\n', 'prior E2f', 'HT xi=2', 'HTHewitt', 'Puck');
for i = 1:size(priors, 1)
    pri = [priors(i,1) priors(i,2) 3 1 1.5 1];
    rand('state', 31); randn('state', 31);
    [pm, ms] = rjmcmc_quick(x, y, pri, prop, N, jumpPr, [priors(i,1) 3 1.5]);
    fprintf('N(%2d,%2d) | %6.1f%%  %6.1f%%  %6.1f%%\n', priors(i,1), priors(i,2), 100*pm);
    % broj skokova (promjena modela) - dijagnostika mijesanja
    fprintf('            broj promjena modela: %d\n', sum(diff(ms) ~= 0));
end
