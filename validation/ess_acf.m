function [ess, rho1] = ess_acf(x)
%ESS_ACF  Efektivna velicina uzorka (initial positive sequence) i lag-1 autokorelacija.
x = x(:) - mean(x);
n = numel(x);
v = (x' * x) / n;
if v == 0
    ess = n; rho1 = 0; return
end
% autokorelacija preko FFT-a
nfft = 2^nextpow2(2 * n);
fx = fft(x, nfft);
acov = real(ifft(fx .* conj(fx)));
acov = acov(1:n) / (n * v);
rho1 = acov(2);
% initial positive sequence
s = 0;
for k = 1:min(n - 1, 1000)
    if acov(k + 1) <= 0
        break
    end
    s = s + acov(k + 1);
end
tau = 1 + 2 * s;
ess = n / max(tau, 1);
end
