function [samples] = sampleHTrjmcmc(x, y, poc,prop,pri, modelKompozita,nSamples)
% poèetne vrijednosti
if nargin<7
    nSamples = 1;
end

% poc=[E2v_0,Em_0,s_0];
c1=poc(1); %odgovara E2v
c2=poc(2); %odgovara Em
c3=poc(3); %odgovara s
% parametri za normal proposal distribuciju
c1Sigma=prop(1); 
c2Sigma=prop(2);
c3Sigma=prop(3);

% % prior istribucija neka bude gamma!
% priorDist=@(x,a,b)gampdf(x,a,1/b);

samples = zeros(nSamples,3);
for si = 1:nSamples
    pc1 = c1 + randn*c1Sigma; % proposal koji genijelac, kako jednostavno !!! vridi samo kada je normalna!!!
    pc2 = c2 + randn*c2Sigma;
    pc3 = c3 + randn*c3Sigma;

    %% dalje nastavi kada budeš ima koncentracije!!! 29.7.22
    pr_b=normpdf(pc1,pri(1) ,pri(2)).*normpdf(c2,pri(3) ,pri(4)).*normpdf(c3,pri(5) ,pri(6)); % ovo priori u obliku normalne distribucija !
    pr_n=normpdf(c1 ,pri(1), pri(2)).*normpdf(c2,pri(3) ,pri(4)).*normpdf(c3,pri(5) ,pri(6));
    
   
        metropolisRatio_1 = exp( likelihoodHTrjmcmc(x,y,pc1,c2,c3,modelKompozita))*pr_b./ ...
       exp( likelihoodHTrjmcmc(x,y,c1,c2,c3,modelKompozita) )/pr_n; % acceptance ratio
    
    if rand < min( 1, metropolisRatio_1 )
        c1 = pc1;
        samples(si,1) = pc1;
    else
        samples(si,1) = c1 ;
    end
    %%
    pr_b=normpdf(c1,pri(1) ,pri(2)).*normpdf(pc2,pri(3) ,pri(4)).*normpdf(c3,pri(5) ,pri(6)); % ovo priori u obliku normalne distribucija !
    pr_n=normpdf(c1 ,pri(1), pri(2)).*normpdf(c2,pri(3) ,pri(4)).*normpdf(c3,pri(5) ,pri(6));

       metropolisRatio_2 = exp( likelihoodHTrjmcmc(x,y,c1,pc2,c3,modelKompozita))*pr_b./ ...
       exp( likelihoodHTrjmcmc(x,y,c1,c2,c3,modelKompozita) )/pr_n; % acceptance ratio

    if rand < min( 1, metropolisRatio_2 )
        c2 = pc2;
        samples(si,2) = c2;
    else
        samples(si,2) = c2;
    end
    
%%
    
    pr_b=normpdf(c1,pri(1) ,pri(2)).*normpdf(c2,pri(3) ,pri(4)).*normpdf(pc3,pri(5) ,pri(6));
    pr_n=normpdf(c1 ,pri(1), pri(2)).*normpdf(c2,pri(3) ,pri(4)).*normpdf(c3,pri(5) ,pri(6));
    

       metropolisRatio_3 = exp( likelihoodHTrjmcmc(x,y,c1,c2,pc3,modelKompozita))*pr_b./ ...
       exp( likelihoodHTrjmcmc(x,y,c1,c2,c3,modelKompozita) )/pr_n; % acceptance ratio
 
    if rand < min( 1, metropolisRatio_3 )
        c3 = pc3;
        samples(si,3) = c3;
    else
        samples(si,3) = c3;
    end
    
end
