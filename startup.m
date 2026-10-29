%STARTUP  Add the folders of this package to the MATLAB path.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'models'));
addpath(fullfile(here, 'mcmc'));
addpath(fullfile(here, 'rjmcmc'));
addpath(fullfile(here, 'validation'));
fprintf('Paths added: models, mcmc, rjmcmc, validation\n');
fprintf('Data folder: %s\n', fullfile(here, 'data'));
