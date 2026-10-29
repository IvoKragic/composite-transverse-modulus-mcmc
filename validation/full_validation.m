%FULL_VALIDATION  Kompletna simulacijska studija za reviziju rada (v0059).
%
%   matlab -batch "run('C:/Users/imari/Documents/StipeRad/tools/full_validation.m')"
%
%   Dijelovi:
%     1) 20 replika po razini suma -> bias, RMSE, pokrivenost 95% intervala,
%        tocnost odabira modela (s prikladnim priorom na sigma)
%     2) ponovljivost odabira modela kroz 5 seedova x 5 priora za E2f
%     3) dijagnostika (stopa prihvata, ESS, lag-1) za tri sinteticka slucaja
%
%   Rezultati: docs/proracuni/validacija.csv i ispis u konzolu.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
cd(root); addpath('models'); addpath('mcmc'); addpath('tools');
outdir = fullfile('docs', 'proracuni');
if ~exist(outdir, 'dir'), mkdir(outdir); end

Vf = linspace(0.3, 0.8, 12);
TRUE = [65 4];                     % E2f, Em
pri  = [55 10  3 0.5  5 5];        % prior: sirok za E2f, pošten za sigma
prop = [1 0.1 0.1];
Nrep = 20; Nmc = 30000; Nrj = 30000; jumpPr = 0.9;

levels = [1.2 ; 2.0 ; 8.0];
truemodel = { 'original', 'Hewitt', 'Hewitt' };

R = {};   % {razina, replika, biasE2f, rmseE2f, covE2f, ...}

fprintf('=== 1) Simulacijska studija: %d replika po razini ===\n', Nrep);
fprintf('%6s | %-22s | %-22s | %-22s | %s\n', ...
    'sigma', 'E2f bias/RMSE/cov', 'Em bias/RMSE/cov', 'sigma bias/RMSE/cov', 'model %');

for L = 1:numel(levels)
    sg = levels(L); md = truemodel{L};
    E2f_est = zeros(Nrep,1); Em_est = zeros(Nrep,1); s_est = zeros(Nrep,1);
    cov = zeros(Nrep,3); okmodel = zeros(Nrep,1);

    for r = 1:Nrep
        seed = 1000*L + r;
        rand('state', seed); randn('state', seed);
        if strcmp(md,'Puck'), E2t = E_2Puck(Vf, TRUE(1), TRUE(2), 0.35);
        else, E2t = E_2HalpinTsai(Vf, TRUE(1), TRUE(2), md); end
        y = E2t + sg*randn(1, numel(Vf));

        rand('state', seed + 7); randn('state', seed + 7);
        ch = rw_mcmc(Vf, y, [55 3 2], prop, pri, md, Nmc);
        s = ch(round(0.1*Nmc)+1:end, :);
        E2f_est(r) = mean(s(:,1)); Em_est(r) = mean(s(:,2)); s_est(r) = mean(s(:,3));
        q = quantile(s, [0.025 0.975]);
        cov(r,1) = q(1,1) <= TRUE(1) && TRUE(1) <= q(2,1);
        cov(r,2) = q(1,2) <= TRUE(2) && TRUE(2) <= q(2,2);
        cov(r,3) = q(1,3) <= sg      && sg      <= q(2,3);

        rand('state', seed + 13); randn('state', seed + 13);
        pm = rjmcmc_quick(Vf, y, pri, prop, Nrj, jumpPr, [55 3 2]);
        [~, im] = max(pm);
        okmodel(r) = im == find(strcmp({'original','Hewitt','Puck'}, md));
    end

    b = [mean(E2f_est)-TRUE(1), mean(Em_est)-TRUE(2), mean(s_est)-sg];
    rmse = [sqrt(mean((E2f_est-TRUE(1)).^2)), sqrt(mean((Em_est-TRUE(2)).^2)), sqrt(mean((s_est-sg).^2))];
    cv = mean(cov, 1);
    fprintf('%6.1f | %+6.1f / %5.1f / %4.0f%% | %+6.2f / %5.2f / %4.0f%% | %+6.2f / %5.2f / %4.0f%% | %3.0f%%\n', ...
        sg, b(1), rmse(1), 100*cv(1), b(2), rmse(2), 100*cv(2), b(3), rmse(3), 100*cv(3), 100*mean(okmodel));
    R(end+1, :) = {sg, 'summary', b(1), rmse(1), cv(1), b(2), rmse(2), cv(2), b(3), rmse(3), cv(3), mean(okmodel)}; %#ok<SAGROW>
end

fprintf('\n  (cov = udio replika gdje je istinita vrijednost unutar 95%% intervala; model %% = udio replika gdje je odabran pravi model)\n');

%% ================= 2) ponovljivost odabira modela ==================
fprintf('\n=== 2) Stabilnost odabira modela (5 seedova x 5 priora), eksperimentalni podaci ===\n');
d = importdata('E2podaci.dat'); x = d(:,1)'; y = d(:,2)';
priors = [73 5 ; 80 5 ; 65 15 ; 60 15 ; 55 15];
races = [0.5 0.05 0.05]; Nrj2 = 60000;
fprintf('%9s | %-24s | skokovi\n', 'prior E2f', 'HT (5 seedova: mean +- sd)');
for i = 1:size(priors,1)
    pp = [priors(i,1) priors(i,2) 3 1 1.5 1];
    ph = zeros(5,1); nj = zeros(5,1);
    for k = 1:5
        rand('state', 100+k); randn('state', 100+k);
        [pm, ms] = rjmcmc_quick(x, y, pp, races, Nrj2, jumpPr, [priors(i,1) 3 1.5]);
        ph(k) = pm(1); nj(k) = sum(diff(ms) ~= 0);
    end
    fprintf('N(%2d,%2d) | %5.1f%% +- %4.1f%%  (min %4.1f, max %4.1f) | %d\n', ...
        priors(i,1), priors(i,2), 100*mean(ph), 100*std(ph), 100*min(ph), 100*max(ph), round(mean(nj)));
end

%% ================= 3) dijagnostika =================================
fprintf('\n=== 3) Dijagnostika lanca (sinteticki slucajevi, N=50000) ===\n');
for L = 1:numel(levels)
    sg = levels(L); md = truemodel{L};
    seed = 500 + L;
    rand('state', seed); randn('state', seed);
    if strcmp(md,'Puck'), E2t = E_2Puck(Vf, TRUE(1), TRUE(2), 0.35);
    else, E2t = E_2HalpinTsai(Vf, TRUE(1), TRUE(2), md); end
    y = E2t + sg*randn(1, numel(Vf));
    rand('state', seed+3); randn('state', seed+3);
    [ch, acc] = rw_mcmc(Vf, y, [55 3 1], prop, pri, md, 50000);
    s = ch(5001:end, :);
    for j = 1:3
        [ess, rho1] = ess_acf(s(:,j));
        fprintf('  sigma=%4.1f par%d: acc=%4.1f%%  ESS=%5.0f  rho1=%.2f\n', sg, j, 100*acc(j), ess, rho1);
    end
end

T = cell2table(R, 'VariableNames', {'sigma','tip','bias_E2f','rmse_E2f','cov_E2f', ...
    'bias_Em','rmse_Em','cov_Em','bias_sigma','rmse_sigma','cov_sigma','odabir_modela'});
writetable(T, fullfile(outdir, 'validacija.csv'));
fprintf('\nSpremljeno: %s\n', fullfile(outdir,'validacija.csv'));
