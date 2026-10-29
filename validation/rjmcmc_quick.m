function [pModel, modelSamples] = rjmcmc_quick(x, y, pri, prop, nSamples, jumpPr, poc)
%RJMCMC_QUICK  Kratki model-jumping lanac (jednake dimenzionalnosti).
%
%   Sva tri modela dijele isti parametarski prostor [E2f Em sigma], pa je
%   prijelaz izmedju modela obican MH korak na diskretnom prostoru modela:
%   prijedlog je uniforman medju preostala dva modela, parametri se ne
%   mijenjaju, prior i Jacobian se pokrate, a prihvat je min(1, L'/L).
%
%   Ulaz:  pri = [mu1 sd1 mu2 sd2 mu3 sd3], prop = [s1 s2 s3],
%          jumpPr = vjerojatnost pokusaja skoka, poc = [E2f Em sigma]
%   Izlaz: pModel = [p1 p2 p3] (udio uzoraka po modelu, nakon 10% burn-in)

models = {'original', 'Hewitt', 'Puck'};
m = 1;
c = poc(:)';
ll = likelihoodHT(x, y, c(1), c(2), c(3), models{m});

    function p = prior(a, b, s)
        p = normpdf(a, pri(1), pri(2)) * normpdf(b, pri(3), pri(4)) * ...
            normpdf(s, pri(5), pri(6));
        if s <= 0, p = 0; end
    end

pr = prior(c(1), c(2), c(3));
modelSamples = zeros(nSamples, 1);

for it = 1:nSamples
    if rand < jumpPr
        % --- skok na drugi model (parametri ostaju) ---
        drugi = setdiff(1:3, m);
        mn = drugi(randi(2));
        lln = likelihoodHT(x, y, c(1), c(2), c(3), models{mn});
        if rand < min(1, exp(lln - ll))
            m = mn; ll = lln;
        end
    else
        % --- MH sweep unutar modela ---
        for j = 1:3
            cn = c; cn(j) = c(j) + prop(j) * randn;
            if cn(3) <= 0, continue; end
            lln = likelihoodHT(x, y, cn(1), cn(2), cn(3), models{m});
            prn = prior(cn(1), cn(2), cn(3));
            if rand < min(1, exp(lln - ll) * prn / pr)
                c = cn; ll = lln; pr = prn;
            end
        end
    end
    modelSamples(it) = m;
end

s = modelSamples(round(0.1 * nSamples) + 1:end);
pModel = [mean(s == 1), mean(s == 2), mean(s == 3)];
end
