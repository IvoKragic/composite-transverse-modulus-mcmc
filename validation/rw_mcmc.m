function [chain, acc] = rw_mcmc(x, y, poc, prop, pri, model, nSamples)
%RW_MCMC  Random-walk Metropolis, jedan parametar po koraku (kao sampleHT).
%
%   Isti algoritam kao mcmc/sampleHT.m, ali:
%     - prior na sigma je odsjecen na sigma > 0 (fizikalno jedino valjano),
%     - vraca broj prihvata po parametru (za dijagnostiku).
%
%   Ulaz:
%     x, y      podaci (Vf, E2)
%     poc       pocetne vrijednosti [E2f Em sigma]
%     prop      standardne devijacije prijedloga [E2f Em sigma]
%     pri       [mu1 sd1 mu2 sd2 mu3 sd3] za E2f, Em, sigma
%     model     'original' | 'Hewitt' | 'Puck'
%     nSamples  broj iteracija
%   Izlaz:
%     chain     [nSamples x 3]
%     acc       [1 x 3] udio prihvata

c1 = poc(1); c2 = poc(2); c3 = poc(3);

% zajednicki prior (nezavisni normalni), sigma odsjecena na > 0
    function p = prior(a, b, c)
        p = normpdf(a, pri(1), pri(2)) * normpdf(b, pri(3), pri(4)) * ...
            normpdf(c, pri(5), pri(6));
        if c <= 0
            p = 0;
        end
    end

chain = zeros(nSamples, 3);
acc = zeros(1, 3);
ll0 = likelihoodHT(x, y, c1, c2, c3, model);
pr0 = prior(c1, c2, c3);

for si = 1:nSamples
    % --- E2f ---
    pc1 = c1 + prop(1) * randn;
    if pc1 > 0
        ll = likelihoodHT(x, y, pc1, c2, c3, model);
        prb = prior(pc1, c2, c3);
        if rand < min(1, exp(ll - ll0) * prb / pr0)
            c1 = pc1; ll0 = ll; pr0 = prb; acc(1) = acc(1) + 1;
        end
    end
    % --- Em ---
    pc2 = c2 + prop(2) * randn;
    if pc2 > 0
        ll = likelihoodHT(x, y, c1, pc2, c3, model);
        prb = prior(c1, pc2, c3);
        if rand < min(1, exp(ll - ll0) * prb / pr0)
            c2 = pc2; ll0 = ll; pr0 = prb; acc(2) = acc(2) + 1;
        end
    end
    % --- sigma ---
    pc3 = c3 + prop(3) * randn;
    if pc3 > 0
        ll = likelihoodHT(x, y, c1, c2, pc3, model);
        prb = prior(c1, c2, pc3);
        if rand < min(1, exp(ll - ll0) * prb / pr0)
            c3 = pc3; ll0 = ll; pr0 = prb; acc(3) = acc(3) + 1;
        end
    end

    chain(si, :) = [c1 c2 c3];
end
acc = acc / nSamples;
end
