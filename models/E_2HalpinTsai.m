function E2 = E_2HalpinTsai(Vv,E2v,Em,model)
% racuna modul elasticnosti u poprecnom smjeru
% Ulazne varijable:
% Vv - volumenski udjel vlakana
% E2v - poprecni modul elasticnosti vlakana
% Em - modul elasticnosti izotropne matrice
% model - model za racunanje parametra ksi ('kvadrat', omjer w_v/t_v, 'Hewitt', 'Giner')
%
if ischar(model)
    if strcmp(model,'original') % kvadratni raspored kruznih vlakana
        ksi_E2 = 2*ones(1,length(Vv));
    elseif strcmp(model,'Hewitt') % Hewitt-Mahlerbe izraz
        ksi_E2 = 2 + 40*Vv.^10;
    elseif strcmp(model,'Giner') %Ginerov izraz
        for n = 1:length(Vv)
            if Vv(n) < 0.3
                ksi_E2(n) = 4.924-35.888*Vv(n)+125.118*Vv(n)^2-145.212*Vv(n)^3;
            else
                ksi_E2(n) = 1.5+5500*Vv(n)^18;
            end
        end
    end
else % pravokutna vlakna dimenzija w_v*t_v u šesterokutnom rasporedu
    ksi_E2 = 2*model;
end
% parametar eta_E2
eta_E2 = (E2v/Em-1)./(E2v/Em+ksi_E2);
% E_2 prema Halpin-Tsai
E2 = Em*(1+ksi_E2.*eta_E2.*Vv)./(1-eta_E2.*Vv);
