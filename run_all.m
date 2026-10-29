function run_all(mode)
%RUN_ALL  Reproduce Table 6 and Supplementary Tables S1-S4 of the manuscript.
%
%   run_all            quick settings (about 2-3 min)
%   run_all('full')    settings used in the paper (about 15-25 min)
%
%   Block 1  Table 6      experimental data: single-stage MCMC and the second (comparison) stage
%   Block 2  Table S1     effect of the prior on the noise standard deviation
%   Block 3  Table S2     replicate study (bias, RMSE, coverage, model identification)
%   Block 4  Table S4     posterior model probability vs prior on E2f

if nargin < 1
    mode = 'quick';
end
here = fileparts(mfilename('fullpath'));
cd(here);
startup;

switch lower(mode)
    case 'full'
        N = 100000;  Nrep = 20;  Nrj = 30000;
    otherwise
        N = 20000;   Nrep = 5;   Nrj = 10000;
        fprintf('QUICK MODE: values differ slightly from the paper (Monte Carlo error).\n');
        fprintf('Run run_all(''full'') for the published settings.\n');
end
fprintf('\n');

prop = [0.5 0.05 0.05];
d = importdata(fullfile(here, 'data', 'E2podaci.dat'));
x = d(:, 1)';  y = d(:, 2)';
burn = round(0.1 * N);

%% ------------------------------------------------ Block 1: Table 6
fprintf('=== Block 1: experimental data (Table 6) ===\n');
pri_rj  = [73 5   3 1    1.5 1  ];
pri_two = [72 3   3.2 0.5  1  0.3];
for k = 1:2
    if k == 1
        pri = pri_rj;  lbl = 'single-stage MCMC (RJMCMC priors)';
    else
        pri = pri_two; lbl = 'second stage (posterior-informed priors)';
    end
    rand('state', 11); randn('state', 11);
    [ch, acc] = rw_mcmc(x, y, [70 3 1], prop, pri, 'original', N);
    s = ch(burn+1:end, :);
    q = quantile(s, [0.025 0.975]);
    fprintf('  %-40s E2f=%6.2f [%6.2f %6.2f]  Em=%5.2f  sigma=%5.2f [%4.2f %4.2f]  acc=%.0f%%\n', ...
        lbl, mean(s(:,1)), q(1,1), q(2,1), mean(s(:,2)), mean(s(:,3)), q(1,3), q(2,3), 100*acc(1));
end
fprintf('  published: RJMCMC 72.16 / 3.21 / 1.97; single stage 72.19 / 3.22 / 1.96\n');
fprintf('\n');

%% ------------------------------------------------ Block 2: Table S1
fprintf('=== Block 2: prior on the noise standard deviation (Table S1) ===\n');
Vf = linspace(0.3, 0.8, 12);
mdl  = {'original', 'Hewitt', 'Hewitt'};
truesig = [1.2 2.0 8.0];
seeds   = [55 1 1];
pri_old = [55 10  3 0.5  0.5 0.3];
pri_new = [55 10  3 0.5  5   5  ];
fprintf('    true sigma | narrow prior | N+(5,5) prior | E2f with N+(5,5)\n');
for i = 1:3
    rand('state', seeds(i)); randn('state', seeds(i));
    ysyn = E_2HalpinTsai(Vf, 65, 4, mdl{i}) + truesig(i)*randn(1, numel(Vf));
    res = zeros(1, 2);  e2f = 0;
    for k = 1:2
        if k == 1, p = pri_old; else, p = pri_new; end
        rand('state', 1000+seeds(i)); randn('state', 1000+seeds(i));
        ch = rw_mcmc(Vf, ysyn, [55 3 0.5], [1 0.1 0.1], p, mdl{i}, N);
        s = ch(burn+1:end, :);
        res(k) = mean(s(:,3));
        if k == 2, e2f = mean(s(:,1)); end
    end
    fprintf('    %10.1f | %12.2f | %13.2f | %6.1f\n', truesig(i), res(1), res(2), e2f);
end
fprintf('  published (narrow / N+(5,5)): 0.95/1.27, 1.19/1.83, 2.50/6.45 GPa\n');
fprintf('\n');

%% ------------------------------------------------ Block 3: Table S2
fprintf('=== Block 3: %d-replicate study (Table S2) ===\n', Nrep);
pri = [55 10 3 0.5 5 5];
fprintf('    true sigma |  E2f bias/RMSE/cov  |  Em bias/RMSE/cov  |  sigma bias/RMSE/cov  | model\n');
for L = 1:3
    E = zeros(Nrep, 3);  cov = zeros(Nrep, 3);  ok = 0;
    for r = 1:Nrep
        seed = 1000*L + r;
        rand('state', seed); randn('state', seed);
        ysyn = E_2HalpinTsai(Vf, 65, 4, mdl{L}) + truesig(L)*randn(1, numel(Vf));
        rand('state', seed+7); randn('state', seed+7);
        ch = rw_mcmc(Vf, ysyn, [55 3 2], [1 0.1 0.1], pri, mdl{L}, N);
        s = ch(burn+1:end, :);
        E(r, :) = mean(s);
        q = quantile(s, [0.025 0.975]);
        cov(r,1) = q(1,1) <= 65 && 65 <= q(2,1);
        cov(r,2) = q(1,2) <= 4  && 4  <= q(2,2);
        cov(r,3) = q(1,3) <= truesig(L) && truesig(L) <= q(2,3);
        rand('state', seed+13); randn('state', seed+13);
        pm = rjmcmc_quick(Vf, ysyn, pri, [1 0.1 0.1], Nrj, 0.9, [55 3 2]);
        [~, im] = max(pm);
        ok = ok + (im == find(strcmp({'original','Hewitt','Puck'}, mdl{L})));
    end
    b  = mean(E) - [65 4 truesig(L)];
    rm = sqrt(mean((E - [65 4 truesig(L)]).^2));
    fprintf('    %10.1f | %+5.1f/%4.1f/%4.0f%% | %+5.2f/%4.2f/%4.0f%% | %+5.2f/%4.2f/%4.0f%% | %3.0f%%\n', ...
        truesig(L), b(1), rm(1), 100*mean(cov(:,1)), b(2), rm(2), 100*mean(cov(:,2)), ...
        b(3), rm(3), 100*mean(cov(:,3)), 100*ok/Nrep);
end
fprintf('  (published with 20 replicates: coverage 80-100%%, correct model 80-95%%)\n');
fprintf('\n');

%% ------------------------------------------------ Block 4: Table S4
fprintf('=== Block 4: model probability vs prior on E2f (Table S4) ===\n');
priors = [73 5; 65 15; 55 15];
for i = 1:size(priors, 1)
    p = [priors(i,1) priors(i,2) 3 1 1.5 1];
    rand('state', 31); randn('state', 31);
    pm = rjmcmc_quick(x, y, p, prop, Nrj, 0.9, [priors(i,1) 3 1.5]);
    fprintf('    N(%2d,%2d) | HT xi=2: %5.1f%%   HT Hewitt: %4.1f%%   Puck: %4.1f%%\n', ...
        priors(i,1), priors(i,2), 100*pm);
end
fprintf('  published: 95.9 +/- 0.6 (N(73,5)), 89.3 (N(65,15)), 80.8 +/- 11.2 (N(55,15))\n');
fprintf('\nDone.\n');
end
