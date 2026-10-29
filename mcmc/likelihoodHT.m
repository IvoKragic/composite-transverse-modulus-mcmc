function loglik = likelihoodHT(Vv,E2_mjereno,E2v,Em,s,model)
n=numel(Vv);
% (x,y,pc1,c2,c3, c4,)
%  E2 = E_2HalpinTsai(Vv,E2v,Em,model)
if strcmp(model,'Puck') % kvadratni raspored kruznih vlakana
    % ovde ubavi pucka!!
    %      E_2Puck(Vv,E2v,Em,NUm)
    
    % Poissonov koeficijent našao u Vlakovim kodovima!! 
            % dilema? dodatan parametar ili uzeti kao gotovu vrijednost??
            % za sada neka bude 0.35
     E2_model=E_2Puck(Vv,E2v,Em,0.35);% ovo mi je model za E2

else 
 E2_model=E_2HalpinTsai(Vv,E2v,Em,model);% ovo mi je model za E2

end

loglik =-n/2*log(2*pi)-n/2*log(s^2)-(E2_mjereno - E2_model)*(E2_mjereno - E2_model)'/(2*s^2); % moje možda se može i bez rvog èlana jer ne utjeèe na oblik nego samo na visinu!!
