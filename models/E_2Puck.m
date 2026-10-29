function E2 = E_2Puck(Vv,E2v,Em,NUm)
% racuna modul elasticnosti u poprecnom smjeru (Puckov mikromehanicki model)
% Ulazne varijable:
% Vv  - volumenski udjel vlakana
% E2v - poprecni modul elasticnosti vlakana
% Em  - modul elasticnosti izotropne matrice
% NUm - Poissonov koeficijent izotropne matrice
%       (ako se ne preda, koristi se 0.35 - vrijednost iz Vlakovih kodova)
%
% Napomena: objedinjene su dvije izvorne verzije funkcije. MCMC pipeline poziva
% s 4 argumenta (NUm=0.35), a RJMCMC pipeline s 3 argumenta. nargin guard daje
% identican rezultat za oba slucaja.
if nargin < 4
    NUm = 0.35;
end
Vm = 1 - Vv;            % volumenski udjel matrice
EmMod = Em/(1-NUm^2);   % modificirani modul elasticnosti matrice
% E_2 prema Pucku
E2 = EmMod*(1+0.85*Vv.^2)./(Vm.^1.25+Vv*EmMod/E2v);
end
