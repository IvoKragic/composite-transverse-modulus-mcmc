%JCM_REVISION_RUNS  Proracuni za odgovor recenzentu (rad v0059, JCM).
%
%   Pokretanje:
%       matlab -batch "run('C:/Users/imari/Documents/StipeRad/tools/jcm_revision_runs.m')"
%
%   Racuna:
%     A) utjecaj priora na sigma u tri sinteticka slucaja (Tab. 1, 3, 4)
%     B) eksperimentalni slucaj s originalnim vs. "informiranim" priorima
%        (provjera dvostrukog koristenja podataka)
%     C) dijagnostika: stope prihvata, ESS, lag-1 autokorelacija, rasipanje
%        izmedju seedova
%     D) nelinearni najmanji kvadrati (NLS) s intervalima pouzdanosti
%
%   Rezultati: docs/proracuni/rezultati.csv (+ ispis u konzolu)

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
cd(root);
addpath('models'); addpath('mcmc'); addpath('tools');
outdir = fullfile('docs', 'proracuni');
if ~exist(outdir, 'dir'), mkdir(outdir); end

fprintf('Radni direktorij: %s\n\n', pwd);
rows = {};   % {grupa, slucaj, prior, seed, par, mean, sd, lo, hi, acc, ess, rho1}

%% ================= A) sinteticki slucajevi ==========================
Vf = linspace(0.3, 0.8, 12);
cases = { ...
    'sigma=1.2 (HT xi=2)',  'original', 1.2, 55 ; ...
    'sigma=2   (HTHewitt)', 'Hewitt',   2.0,  1 ; ...
    'sigma=8   (HTHewitt)', 'Hewitt',   8.0,  1 };
pri_old = [55 10  3 0.5  0.5 0.3];
pri_new = [55 10  3 0.5  5   5  ];    % sigma ~ N(5,5), odsjeceno na > 0
poc = [55 3 0.5]; prop = [1 0.1 0.1]; N = 50000; burn = 0.1 * N;

fprintf('=== A) sinteticki slucajevi: prior na sigma ===\n');
for i = 1:size(cases, 1)
    lbl = cases{i,1}; model = cases{i,2}; sg = cases{i,3}; seed = cases{i,4};

    rand('state', seed); randn('state', seed);
    if strcmp(model, 'Puck')
        E2true = E_2Puck(Vf, 65, 4, 0.35);
    else
        E2true = E_2HalpinTsai(Vf, 65, 4, model);
    end
    y = E2true + sg * randn(1, numel(Vf));

    priornames = {'N(0.5,0.3) - originalni', 'N(5,5) odsjecen na >0 - predlozeni'};
    priors = {pri_old, pri_new};
    for k = 1:2
        rand('state', 1000 + seed); randn('state', 1000 + seed);
        [ch, acc] = rw_mcmc(Vf, y, poc, prop, priors{k}, model, N);
        s = ch(burn+1:end, :);
        fprintf('\n  %s | %s\n', lbl, priornames{k});
        for j = 1:3
            [ess, rho1] = ess_acf(s(:, j));
            q = quantile(s(:, j), [0.025 0.975]);
            fprintf('    par%d  mean=%8.3f  sd=%6.3f  95%%CI=[%7.3f %7.3f]  acc=%5.1f%%  ESS=%6.0f  rho1=%.2f\n', ...
                j, mean(s(:, j)), std(s(:, j)), q(1), q(2), 100*acc(j), ess, rho1);
            rows(end+1, :) = {lbl, priornames{k}, 1000+seed, j, mean(s(:,j)), std(s(:,j)), q(1), q(2), acc(j), ess, rho1}; %#ok<SAGROW>
        end
    end
end

%% ================= B) eksperimentalni podaci =======================
d = importdata('E2podaci.dat');
x = d(:, 1)'; y = d(:, 2)';
fprintf('\n\n=== B) eksperimentalni podaci (n=%d) ===\n', numel(x));

pri_orig = [73 5   3 1    1.5 1  ];   % priore iz RJMCMC analize (kako pise u radu)
pri_inf  = [72 3   3.2 0.5  1  0.3];   % priore "informirane" RJMCMC-om (drugi stupanj)
prop_e = [0.5 0.05 0.05]; N2 = 100000; burn2 = 0.1 * N2;
seeds = [11 12 13];

for k = 1:2
    if k == 1, pri = pri_orig; lbl = 'originalni N(73,5)/N(3,1)/N(1.5,1)';
    else,      pri = pri_inf;  lbl = 'informirani N(72,3)/N(3.2,0.5)/N(1,0.3)';
    end
    M = []; A = [];
    for s = seeds
        rand('state', s); randn('state', s);
        [ch, acc] = rw_mcmc(x, y, [70 3 1], prop_e, pri, 'original', N2);
        M(end+1, :) = mean(ch(burn2+1:end, :)); %#ok<SAGROW>
        if s == seeds(1)
            smp = ch(burn2+1:end, :);
            for j = 1:3
                [ess, rho1] = ess_acf(smp(:, j));
                q = quantile(smp(:, j), [0.025 0.975]);
                fprintf('\n  %s | par%d  mean=%7.3f sd=%5.3f 95%%CI=[%6.3f %6.3f] acc=%5.1f%% ESS=%5.0f rho1=%.2f\n', ...
                    lbl, j, mean(smp(:,j)), std(smp(:,j)), q(1), q(2), 100*acc(j), ess, rho1);
                rows(end+1, :) = {['eksperiment: ' lbl], 'all seeds', s, j, mean(smp(:,j)), std(smp(:,j)), q(1), q(2), acc(j), ess, rho1}; %#ok<SAGROW>
            end
        end
    end
    fprintf('  %s -> srednja vrijednost po seedovima: E2f=%.3f (sd medju seedovima %.3f), Em=%.3f (%.3f), sigma=%.3f (%.3f)\n', ...
        lbl, mean(M(:,1)), std(M(:,1)), mean(M(:,2)), std(M(:,2)), mean(M(:,3)), std(M(:,3)));
end

%% ================= C) NLS usporedba ================================
fprintf('\n\n=== C) nelinearni najmanji kvadrati (NLS) na eksperimentalnim podacima ===\n');
models = {'original', 'Hewitt', 'Puck'};
n = numel(x);
for i = 1:numel(models)
    md = models{i};
    f = @(b, xx) local_model(xx, b(1), b(2), md);
    b0 = [70 3.5];
    opts = statset('Display', 'off', 'MaxIter', 2000, 'TolFun', 1e-10, 'TolX', 1e-10);
    [b, R, ~, covb] = nlinfit(x, y, f, b0, opts);
    se = sqrt(diag(covb));
    rss = sum(R.^2);
    sMLE = sqrt(rss / (n - 2));
    bic = n * log(rss / n) + 2 * log(n);
    fprintf('  %-8s  E2f=%7.2f (95%%CI %6.2f-%6.2f)  Em=%5.2f (95%%CI %5.2f-%5.2f)  sigma_MLE=%.2f  RSS=%.1f  BIC=%.1f\n', ...
        md, b(1), b(1)-1.96*se(1), b(1)+1.96*se(1), b(2), b(2)-1.96*se(2), b(2)+1.96*se(2), sMLE, rss, bic);
end

%% ================= spremi ==========================================
T = cell2table(rows, 'VariableNames', ...
    {'grupa','prior','seed','parametar','mean','sd','ci_lo','ci_hi','acc','ess','rho1'});
writetable(T, fullfile(outdir, 'rezultati.csv'));
fprintf('\n\nSpremljeno: %s\n', fullfile(outdir, 'rezultati.csv'));

function E2 = local_model(Vv, E2v, Em, md)
if strcmp(md, 'Puck')
    E2 = E_2Puck(Vv, E2v, Em, 0.35);
else
    E2 = E_2HalpinTsai(Vv, E2v, Em, md);
end
end
