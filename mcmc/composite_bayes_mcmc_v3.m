clear all ;clc ;close all

podaciExp =importdata('E2podaci.dat');

podaci='vlak';% vlak petra petra2 simulirano uvezeno 
podatakaUK=12;
zanimaMeModel='original' ;% original, Hewitt ; Puck, Giner njega sam izbacio naknadno
crtaj='pojedinacno' %pojedinacno sve
noise=1.2;

switch podaci
    case 'vlak'
        xpodaci=podaciExp(:,1)';
        ypodaci=podaciExp(:,2)';
    case 'petra'
        dataEXP=load('petraExp.txt');
        xpodaci=dataEXP(:,1);
        ypodaci=dataEXP(:,2);
    case 'petra2'
        dataEXP=load('petraExp2.txt');
        xpodaci=dataEXP(:,1);
        ypodaci=dataEXP(:,2);
    case 'simulirano'
        E2v_provjera=65;
        Em_provjera=4;
        Vv_provjera=linspace(0.3,0.8,podatakaUK);
        if strcmp(zanimaMeModel,'Puck')
            %             E2 = E_2Puck(Vv,E2v,Em,NUm)
            E2_provjera=E_2Puck(Vv_provjera,E2v_provjera,Em_provjera,0.35); %samo za Pucka
            zaprikazOriginala=E2_provjera;
            seed=50; rand( 'state' , seed ); randn('state',seed );
            E2_provjera=E2_provjera+normrnd(0,noise,size(E2_provjera));
            
        else
            E2_provjera=E_2HalpinTsai(Vv_provjera,E2v_provjera,Em_provjera,zanimaMeModel); %kvadrat Hewitt Giner
            zaprikazOriginala=E2_provjera;
            seed=55; rand( 'state' , seed ); randn('state',seed );
            E2_provjera=E2_provjera+normrnd(0,noise,size(E2_provjera));
            
        end
        
        xpodaci=Vv_provjera;
        ypodaci=E2_provjera;
    case 'uvezeno'
        brFilea='10' % file 10 je za primjer 1
        nodes= xlsread(['podaciZaRJ_',brFilea,'.xlsx']);
        xpodaci=nodes(:,1)';
        ypodaci=nodes(:,2)';
        if strcmp(brFilea,'10')
            %         ovo u podaci samo za RJ10!
            Vv_provjera=linspace(0.0,1,100);
            E2v_provjera=65;
            Em_provjera=4;
            zaprikazOriginala_RJ=E_2HalpinTsai(Vv_provjera,E2v_provjera,Em_provjera,'original');
        elseif strcmp(brFilea,'12')
            %         ovo u podaci samo za RJ12!
            Vv_provjera=linspace(0.0,1,100);
            E2v_provjera=65;
            Em_provjera=4;
            zaprikazOriginala_RJ=E_2HalpinTsai(Vv_provjera,E2v_provjera,Em_provjera,'Hewitt');
        elseif strcmp(brFilea,'13')
            %         ovo u podaci samo za RJ12!
            Vv_provjera=linspace(0.0,1,100);
            E2v_provjera=65;
            Em_provjera=4;
            zaprikazOriginala_RJ=E_2HalpinTsai(Vv_provjera,E2v_provjera,Em_provjera,'Hewitt');
        end
%         slejdeæi case možda izbaciti   
end

%%

njihovFitVv=0.2:0.01:1;
njihovFit=E_2HalpinTsai(njihovFitVv,73.1,3.45,'original');% ovo mi je model za E2

% % epoksi - izotropna matrica
pois = 0.35; % Poissonov koeficijent našao u Vlakovim kodovima!!
% dilema? dodatan parametar ili uzeti kao gotovu vrijednost???
%%
if ~strcmp(podaci,'simulirano')
    %     E2_model_Puck=E_2Puck(njihovFitVv,73.1,3.45,pois);% ovo mi je model za E2 usoreðujem Pucka sa njihovim moelima
    %     FS=18;
    %     Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['njihov fit' ],...
    %         'NumberTitle','off','Position',[0.1 0.1 0.8 0.72]);
    %     ax0= axes('Parent',Fig_Optimist,'Position',[0.1 0.15 0.83 0.75],...
    %         'FontName','Times New Roman','FontSize',FS,...
    %         'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5);hold on,
    %     plot(ax0,xpodaci,ypodaci,'sr','markersize',10, 'linewidth',2)
    %     hold on
    %     plot(ax0,njihovFitVv,njihovFit, 'k','linewidth',3)
    %
    %     xlabel(ax0,['Fiber volume','\it V_f '],'FontSize',FS);ylabel(['Elastic modulus','{\it E_2}', ' [GPa]'],'FontSize',FS)
    %
    %     grid on
    
end
%% za raèun prikaz
E2v_linspace=linspace(1e-1,100,1e3);
Em_linspace=linspace(1e-1,10,1e3);
s_raspon=linspace(1e-3,2,1e3);

%% MCMC
N_mcmc=1e5;% broj iteracija
% postavke za prvu predloženu vrijednost
seed=50; rand( 'state' , seed ); randn('state',seed );
% poèetne vrijednosti
E2v_0 = unifrnd( min(E2v_linspace) , max(E2v_linspace));
Em_0= unifrnd( min(Em_linspace) , max(Em_linspace));
s_0 = unifrnd( min(s_raspon) , max(s_raspon));

%% MCMC
burnIn=0.1*N_mcmc;

poc=[E2v_0,Em_0,s_0];


if strcmp(podaci, 'vlak')
    % za NORMAL!!!!!!
    a_1= 72; b_1=3 ;
    a_2= 3.2; b_2=0.5 ;
    a_3= 1; b_3=0.3 ;
    
    %parametri za Normal proposal distribuciju
    proposal_E2v=0.5;
    proposal_Em=0.05;
    proposal_s=0.05;
    
else
    %     % za NORMAL!!!!!!
    a_1= 55; b_1=10 ; % kada je mali šum prior oupæe nije bitan!
    a_2= 3; b_2=0.5 ;
    a_3= 0.5; b_3=0.3 ;
    %parametri za Normal proposal distribuciju
    proposal_E2v=1;
    proposal_Em=0.1;
    proposal_s=0.1;
    
end

prop=[proposal_E2v,proposal_Em,proposal_s];
pri=[a_1 b_1 a_2  b_2 a_3 b_3 ];
figure
subplot(3,1,1)
plot(E2v_linspace, normpdf(E2v_linspace,a_1 ,b_1))
title('E_{2f}')
subplot(3,1,2)
plot(Em_linspace, normpdf(Em_linspace,a_2 ,b_2))
title('E_m')
subplot(3,1,3)
plot(s_raspon, normpdf(s_raspon,a_3 ,b_3))
title('\sigma')

switch  crtaj
    case 'pojedinacno'
        FS=18;
        disp(['gledam samo-',zanimaMeModel])
        %============MCMC====================================
        % sampleHT(x, y, poc,prop,nSamples,pri, model)
        tic
        samplesTot = sampleHT(xpodaci,ypodaci,poc,prop ,N_mcmc,pri,zanimaMeModel);
        
        srednjVrijednost= mean(samplesTot(burnIn:end,:));
        toc
        stipePlotHT(N_mcmc,50, samplesTot, srednjVrijednost,burnIn, zanimaMeModel)
        
        %  E2 = E_2HalpinTsai(Vv,E2v,Em,model)
        Vv_model=0.2:0.01:0.85;
        E2v_Bayes=srednjVrijednost(1)
        Em_Bayes=srednjVrijednost(2)
        s_Bayes=srednjVrijednost(3)
        
        if strcmp(zanimaMeModel,'Puck') % kvadratni raspored kruznih vlakana
            % Poissonov koeficijent našao u Vlakovim kodovima!!
            % dilema? dodatan parametar ili uzeti kao gotovu vrijednost??
            E2_model=E_2Puck(Vv_model,E2v_Bayes,Em_Bayes, pois);% ovo mi je model za E2
        else
            E2_model=E_2HalpinTsai(Vv_model,E2v_Bayes,Em_Bayes,zanimaMeModel);% ovo mi je model za E2
        end
        
        koliko=3;
        E2_model_g=E2_model+koliko*srednjVrijednost(3);
        E2_model_d=E2_model-koliko*srednjVrijednost(3);
        
        % omjer=50;
        % eta=(omjer-1)/(omjer+2);
        % em=16*(1-eta*0.6)/(1+2*eta*0.6)
        % e2v=omjer*em
        
        %% ============MATLABOV FIT====================================
        [xData, yData] = prepareCurveData( xpodaci, ypodaci );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [10 10];
        if strcmp(zanimaMeModel,'original') % kvadratni raspored kruznih vlakana
            % Fit: 'original'.
            ksi_E2='2';
            eta_E2=[ '(E2v/Em-1)./(E2v/Em+',ksi_E2,')'];
            ft_kvadrat= fittype( ['Em*(1+2.*',eta_E2,'.*Vv)./(1-',eta_E2,'.*Vv);'], 'independent', 'Vv', 'dependent', 'y' );
            % Fit model kvadratta.
            [fitresult, gof] = fit( xData, yData, ft_kvadrat, opts );
            
            x=[0.2:0.01:0.85];
            YMatlabFit = feval(fitresult, x);
            
            E2v_MatlabFit=fitresult.E2v;
            Em_MatlabFit=fitresult.Em;
            % ci = confint(fitresult,0.95)
            
        elseif strcmp(zanimaMeModel,'Hewitt')
            % Fit: 'Hewit'.
            ksi_E2 = '2 + 40*Vv.^10';
            eta_E2=[ '(E2v/Em-1)./(E2v/Em+',ksi_E2,')'];
            ft_hewit= fittype( ['Em*(1+',ksi_E2,'.*',eta_E2,'.*Vv)./(1-',eta_E2,'.*Vv);'], 'independent', 'Vv', 'dependent', 'y' );
            % Fit model Hewit.
            [fitresult, gof] = fit( xData, yData, ft_hewit, opts );
            x=[0.2:0.01:0.85];
            YMatlabFit = feval(fitresult, x);
            E2v_MatlabFit=fitresult.E2v;
            Em_MatlabFit=fitresult.Em;
            
        elseif strcmp(zanimaMeModel,'Giner')
            % Fit: 'Giner'
            ksi_E2_1 = '(4.924-35.888*Vv+125.118*Vv.^2-145.212*Vv.^3)';
            ksi_E2_2= '(1.5+5500*Vv.^18)';
            eta_E2_1 = ['(E2v/Em-1)./(E2v/Em+',ksi_E2_1,')'];
            eta_E2_2 = ['(E2v/Em-1)./(E2v/Em+',ksi_E2_2,')'];
            funkcija=['(Em.*(1+',ksi_E2_1,'.*',eta_E2_1,'.*Vv)./(1-',eta_E2_1,'.*Vv)).*(Vv < 0.3)+(Em*(1+',ksi_E2_2,'.*',eta_E2_2,'.*Vv)./(1-',eta_E2_2,'.*Vv)).*(Vv >= 0.3);'];
            ft_Giner= fittype(funkcija, 'independent', 'Vv', 'dependent', 'y' );
            
            [fitresult, gof] = fit( xData, yData, ft_Giner, opts );
            
            x=[0.2:0.01:0.85];
            YMatlabFit = feval(fitresult, x);
            
            E2v_MatlabFit=fitresult.E2v;
            Em_MatlabFit=fitresult.Em;
            
        elseif strcmp(zanimaMeModel,'Puck')
            % Fit: 'Puck'
            Vm = '1 - Vv'; % volumenski udjel matrice
            % modificirani modul elasti?nosti matrice
            EmMod = ['Em/(1-',num2str(pois),'^2)'];
            % E_2 prema Pucku
            % E2 = EmMod*(1+0.85*Vv.^2)./(Vm.^1.25+Vv*EmMod/E2v);
            
            ft_Puck= fittype( [EmMod,'*(1+0.85*Vv.^2)./(',Vm,'.^1.25+Vv*',EmMod,'./E2v)';], 'independent', 'Vv', 'dependent', 'y' );
            % Fit model Hewit.
            [fitresult, gof] = fit( xData, yData, ft_Puck, opts );
            x=[0.2:0.01:0.85];
            YMatlabFit = feval(fitresult, x);
            
            E2v_MatlabFit=fitresult.E2v;
            Em_MatlabFit=fitresult.Em;
            
        end
        
        %% ============PLOTANJE====================================
        
        Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['fitano pomoæu Bayesa' ],...
            'NumberTitle','off','Position',[0.1 0.1 0.8 0.72]);
        ax = axes('Parent',Fig_Optimist,'Position',[0.1 0.15 0.83 0.75],...
            'FontName','Times New Roman','FontSize',FS,...
            'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5);hold on,
        x2 = [Vv_model, fliplr(Vv_model)];
        inBetween = [E2_model_g, fliplr(E2_model_d)];
        P3=fill(ax, x2,inBetween,[0.5 0.5 0.5],'EdgeColor' , 'none','FaceAlpha',0.3);
        P1=plot(ax,xpodaci,ypodaci, 'sr','markersize',10, 'linewidth',1);
        hold on
        P2=plot(ax,Vv_model,E2_model,'-','color',[0 0 1], 'linewidth',4);
        hold on
        P6=plot(ax,x, YMatlabFit, ':', 'color',[0 1 0],'linewidth',3);
        %%crta targer
        if strcmp(podaci, 'simulirano')
            P5=plot(ax,xpodaci, zaprikazOriginala, '--', 'color',[0.9 0.1 0.1],'linewidth',2);
            legend(ax,[ P2,P3,P5, P1],  'Bayesian fit', [char(177), num2str(koliko),'\sigma_{LH}'],'Target','Data','location','best',  'FontSize',FS)
        elseif strcmp(podaci, 'uvezeno')
            P5=plot(ax,Vv_provjera, zaprikazOriginala_RJ, '--', 'color',[0.9 0.1 0.1],'linewidth',2);
            legend(ax,[ P2,P3,P5, P1,P6],  'Bayesian fit', [char(177), num2str(koliko),'\sigma_{LH}'],'Target','Data','Matlab fit','location','best',  'FontSize',FS)
    
        else
            P5=plot(ax,njihovFitVv,njihovFit, '--', 'color',[0.9 0.1 0.1],'linewidth',2);
            legend(ax,[ P2,P3,P5, P1],  'HT _{\xi=2,  MCMC}', [char(177), num2str(koliko),'\sigma_{LH}'],'HT _{\xi=2, [4]}','Experiment [4]','location','best',  'FontSize',FS)
        end
        legend boxoff
        %         xlabel(' Fiber volume V_f ','FontSize',FS);ylabel(' Elastic modulus E_{f2}  [GPa]','FontSize',FS)
        xlabel(['Fiber volume','\it V_f '],'FontSize',20);ylabel(['Elastic modulus','{\it E_2}', ' [GPa]'],'FontSize',20)
        xlim(ax,[0.35, 0.81])
        %         ylim(ax,[8 23])
        grid on
        % ============MATLABOV FIT PLOT====================================
        Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['Bayes vs matlab fit funkcija' ],...
            'NumberTitle','off','Position',[0.1 0.1 0.8 0.72]);
        ax = axes('Parent',Fig_Optimist,'Position',[0.1 0.15 0.83 0.75],...
            'FontName','Times New Roman','FontSize',FS,...
            'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5);hold on,
        P6=plot(ax,x, YMatlabFit, '-', 'color',[0.7 0.7 0.7],'linewidth',5);
        hold on
        P1=plot(ax,xpodaci,ypodaci, 'sr','markersize',10, 'linewidth',2);
        hold on
        P2=plot(ax,Vv_model,E2_model, '-','color',[0 0 1], 'linewidth',2);
        hold on
        if strcmp(podaci, 'simulirano')
            P5=plot(ax,xpodaci, zaprikazOriginala, '--', 'color',[0.9 0.1 0.1],'linewidth',2);
            legend(ax,[ P6,P5, P2,P1 ], 'Matlab fit', 'Target','Bayesian fit','Data','location','southeast',  'FontSize',FS)
        elseif strcmp(podaci, 'uvezeno')
            P5=plot(ax,Vv_provjera, zaprikazOriginala_RJ, '--', 'color',[0.9 0.1 0.1],'linewidth',2);
            legend(ax,[ P6,P5, P2,P1 ], 'Matlab fit', 'Target','Bayesian fit','Data','location','southeast',  'FontSize',FS)
        else
            
            P5=plot(ax,njihovFitVv,njihovFit, '-.', 'color',[0.9 0.1 0.1],'linewidth',1);
            legend(ax,[ P6,P5, P2,P1 ], 'Matlab fit', 'Halpin-Tsai','Bayesian fit','Data','location','southeast',  'FontSize',FS)
        end
        
        legend boxoff
        xlabel(['Fiber volume','\it V_f '],'FontSize',20);ylabel(['Elastic modulus','{\it E_2}', ' [GPa]'],'FontSize',20)
        text(ax, 0.36,24,['Matlab fit: E_{2f}=', num2str(E2v_MatlabFit),' E_m=', num2str(Em_MatlabFit)])
        text(ax, 0.36,22,['Bayes fit: E_{2f}=', num2str(E2v_Bayes),' E_m=', num2str(Em_Bayes)])
        
        if strcmp(podaci, 'simulirano')
            text(ax, 0.36,20,['Terget : E_{2f}=',num2str(E2v_provjera),' E_m=',num2str(Em_provjera)])
        elseif strcmp(podaci, 'uvezeno')
            text(ax, 0.36,20,'Terget : E_{2f}=65 E_m=4')
        else
            text(ax, 0.36,20,'Originalno_{clanak} : E_{2f}=73.1 E_m=3.45')
        end
        text(ax, 0.36,18,['data ', num2str(numel(ypodaci))])
        grid on
        xlim([0.35 0.81])
    case 'sve'
        
        samplesTotKvadrat= sampleHT(xpodaci,ypodaci,poc,prop ,N_mcmc,pri,'original');
        samplesTotHewit= sampleHT(xpodaci,ypodaci,poc,prop ,N_mcmc,pri,'Hewitt');
        %         samplesTotGiner= sampleHT(xpodaci,ypodaci,poc,prop ,N_mcmc,pri,'Giner');
        samplesTotPuck= sampleHT(xpodaci,ypodaci,poc,prop ,N_mcmc,pri,'Puck');
        
        srednjVrijednostKvadrat= mean(samplesTotKvadrat(burnIn:end,:));
        srednjVrijednostHewitt= mean(samplesTotHewit(burnIn:end,:));
        %         srednjVrijednostGiner= mean(samplesTotGiner(burnIn:end,:));
        srednjVrijednostPuck= mean(samplesTotPuck(burnIn:end,:));
        
        stipePlotHT(N_mcmc,50, samplesTotKvadrat, srednjVrijednostKvadrat,burnIn,'original' )
        stipePlotHT(N_mcmc,50, samplesTotHewit, srednjVrijednostHewitt,burnIn,'Hewitt')
        %         stipePlotHT(N_mcmc,50, samplesTotGiner, srednjVrijednostGiner,burnIn,'Giner')
        stipePlotHT(N_mcmc,50, samplesTotPuck, srednjVrijednostPuck,burnIn,'Puck')
        
        %%
        %  E2 = E_2HalpinTsai(Vv,E2v,Em,model)
        Vv_model=0.1:0.01:0.8;
        E2v_Bayes_kvadrat=srednjVrijednostKvadrat(1);
        Em_Bayes_kvadrat=srednjVrijednostKvadrat(2);
        s_Bayes_kvadrat=srednjVrijednostKvadrat(3);
        E2_model_kvadrat=E_2HalpinTsai(Vv_model,E2v_Bayes_kvadrat,Em_Bayes_kvadrat,'original');% ovo mi je model za E2
        rezultat_kvadrat=[E2v_Bayes_kvadrat Em_Bayes_kvadrat s_Bayes_kvadrat]
        
        E2v_Bayes_Hewitt=srednjVrijednostHewitt(1);
        Em_Bayes_Hewitt=srednjVrijednostHewitt(2);
        s_Bayes_Hewitt=srednjVrijednostHewitt(3);
        E2_model_Hewitt=E_2HalpinTsai(Vv_model,E2v_Bayes_Hewitt,Em_Bayes_Hewitt,'Hewitt');% ovo mi je model za E2
        rezultat_Hewit=[E2v_Bayes_Hewitt Em_Bayes_Hewitt s_Bayes_Hewitt]
        
        %         E2v_Bayes_Giner=srednjVrijednostGiner(1);
        %         Em_Bayes_Giner=srednjVrijednostGiner(2);
        %         s_Bayes_Giner=srednjVrijednostGiner(3);
        %         E2_model_Giner=E_2HalpinTsai(Vv_model,E2v_Bayes_Giner,Em_Bayes_Giner,'Giner');% ovo mi je model za E2
        %         rezultat_Giner=[E2v_Bayes_Giner Em_Bayes_Giner s_Bayes_Giner]
        
        E2v_Bayes_Puck=srednjVrijednostPuck(1);
        Em_Bayes_Puck=srednjVrijednostPuck(2);
        s_Bayes_Puck=srednjVrijednostPuck(3);
        E2_model_Puck=E_2Puck(Vv_model,E2v_Bayes_Puck,Em_Bayes_Puck,pois);% ovo mi je model za E2
        rezultat_Puck=[E2v_Bayes_Puck Em_Bayes_Puck s_Bayes_Puck]
        
        FS=20;
        Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['fitano pomoæu Bayesa' ],...
            'NumberTitle','off','Position',[0.1 0.1 0.8 0.72]);
        ax1= axes('Parent',Fig_Optimist,'Position',[0.1 0.2 0.83 0.75],...
            'FontName','Times New Roman','FontSize',FS,...
            'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5);hold on,
        P1=plot(ax1,xpodaci,ypodaci, 'sr','markersize',10, 'linewidth',2);
        hold on
        P2=plot(ax1,Vv_model,E2_model_kvadrat, 'k', 'linewidth',2);
        hold on
        P3=plot(ax1,Vv_model,E2_model_Hewitt, '--b', 'linewidth',2);
        %         hold on
        %         P4=plot(ax1,Vv_model,E2_model_Giner, 'c', 'linewidth',2);
        hold on
        P6=plot(ax1,Vv_model,E2_model_Puck, '-.m', 'linewidth',2);
        
        if strcmp(podaci, 'simulirano')
            P5=plot(ax1,xpodaci, zaprikazOriginala, ':', 'color',[0.5 0.5 0.5],'linewidth',2);
            legend(ax1,[ P2 P3 P5 P6 P1],  'HT_{\it \xi=2}','HT_{Hewitt}','Target','Puck','Data', 'location','best','FontSize',FS)
        elseif strcmp(podaci, 'uvezeno')
            P5=plot(ax1,Vv_provjera, zaprikazOriginala_RJ, ':', 'color',[0.5 0.5 0.5],'linewidth',2);
            legend(ax1,[ P2 P3 P5 P6 P1],  'HT_{\it \xi=2}','HT_{Hewitt}','Target','Puck','Data', 'location','best','FontSize',FS)
        else
            P5=plot(ax1,njihovFitVv,njihovFit, ':', 'color',[0.5 0.5 0.5],'linewidth',2);
            legend(ax1,[ P2 P3 P5 P6 P1],  'HT_{\it \xi=2}','HT_{Hewitt}','HT_{[4]}','Puck','Data', 'location','best','FontSize',FS)
        end
        legend boxoff
        xlabel(['Fiber volume','\it V_f '],'FontSize',20);ylabel(['Elastic modulus','{\it E_2}', '[GPa]'],'FontSize',20)
        grid on
        text(ax1, 0.36,28,['HT-Original fit: E_{2v}=', num2str(E2v_Bayes_kvadrat),' E_m=', num2str(Em_Bayes_kvadrat)])
        text(ax1, 0.36,26,['Hewitt fit: E_{2f}=', num2str(E2v_Bayes_Hewitt),' E_m=', num2str(Em_Bayes_Hewitt)])
        %         text(ax1, 0.1,24,['Giner fit: E_{2f}=', num2str(E2v_Bayes_Giner),' E_m=', num2str(E2v_Bayes_Giner)])
        text(ax1, 0.36,24,['Puck fit: E_{2f}=', num2str(E2v_Bayes_Puck),' E_m=', num2str(Em_Bayes_Puck)])
        
        if strcmp(podaci, 'simulirano')
            text(ax1, 0.36,22,['Terget : E_{2f}=',num2str(E2v_provjera),' E_m=',num2str(Em_provjera)])
        elseif strcmp(podaci, 'uvezeno')
            text(ax1, 0.36,22,'Terget : E_{2f}=65 E_m=4')
        else
            text(ax1, 0.36,22,'Originalno_{clanak} : E_{2f}=73.1 E_m=3.45')
        end
        
        xlim(ax1,[0.35,0.81])
        %% ide fitanje po matlabu
        [xData, yData] = prepareCurveData( xpodaci, ypodaci );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [5 0.1];
        %% Fit: 'original'.
        ksi_E2='2';
        eta_E2=[ '(E2v/Em-1)./(E2v/Em+',ksi_E2,')'];
        ft_kvadrat= fittype( ['Em*(1+2.*',eta_E2,'.*Vv)./(1-',eta_E2,'.*Vv);'], 'independent', 'Vv', 'dependent', 'y' );
        % Fit model kvadratta.
        [fitresultKvadrat, gofKvadrat] = fit( xData, yData, ft_kvadrat, opts );
        
        %% Fit: 'Hewit'.
        ksi_E2 = '2 + 40*Vv.^10';
        eta_E2=[ '(E2v/Em-1)./(E2v/Em+',ksi_E2,')'];
        ft_hewit= fittype( ['Em*(1+',ksi_E2,'.*',eta_E2,'.*Vv)./(1-',eta_E2,'.*Vv);'], 'independent', 'Vv', 'dependent', 'y' );
        % Fit model Hewit.
        [fitresultHewit, gofHewit] = fit( xData, yData, ft_hewit, opts );
        %% Fit: 'Giner'
        %         ksi_E2_1 = '(4.924-35.888*Vv+125.118*Vv.^2-145.212*Vv.^3)';
        %         ksi_E2_2= '(1.5+5500*Vv.^18)';
        %         eta_E2_1 = ['(E2v/Em-1)./(E2v/Em+',ksi_E2_1,')'];
        %         eta_E2_2 = ['(E2v/Em-1)./(E2v/Em+',ksi_E2_2,')'];
        %         funkcija=['(Em.*(1+',ksi_E2_1,'.*',eta_E2_1,'.*Vv)./(1-',eta_E2_1,'.*Vv)).*(Vv < 0.3)+(Em*(1+',ksi_E2_2,'.*',eta_E2_2,'.*Vv)./(1-',eta_E2_2,'.*Vv)).*(Vv >= 0.3);'];
        %         ft_Giner= fittype(funkcija, 'independent', 'Vv', 'dependent', 'y' );
        %
        %         [fitresultGiner, gofGiner] = fit( xData, yData, ft_Giner, opts );
        
        %% Fit: 'Puck'
        Vm = '1 - Vv'; % volumenski udjel matrice
        % modificirani modul elasti?nosti matrice
        EmMod = ['Em/(1-',num2str(pois),'^2)'];
        % E_2 prema Pucku
        % E2 = EmMod*(1+0.85*Vv.^2)./(Vm.^1.25+Vv*EmMod/E2v);
        
        ft_Puck= fittype( [EmMod,'*(1+0.85*Vv.^2)./(',Vm,'.^1.25+Vv*',EmMod,'./E2v)';], 'independent', 'Vv', 'dependent', 'y' );
        % Fit model Hewit.
        [fitresultPuck, gofPuck] = fit( xData, yData, ft_Puck, opts );
        
        %% plotanje
        x=[0.3:0.01:0.8];
        YKvadrat = feval(fitresultKvadrat, x);
        YHewit = feval(fitresultHewit, x);
        %         YGiner = feval(fitresultGiner, x);
        YPuck = feval(fitresultPuck, x);
        FS=20;
        Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['fitano pomoæu Matlab FIT funkcije' ],...
            'NumberTitle','off','Position',[0.1 0.1 0.8 0.72]);
        ax2 = axes('Parent',Fig_Optimist,'Position',[0.1 0.15 0.83 0.75],...
            'FontName','Times New Roman','FontSize',FS,...
            'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5);hold on,
        P1=plot(ax2,xpodaci,ypodaci, 'sr','markersize',10, 'linewidth',2);
        hold on
        P2=plot(ax2,x,YKvadrat, 'k', 'linewidth',2);
        hold on
        P3=plot(ax2,x,YHewit, '--b', 'linewidth',2);
        %         hold on
        %         P4=plot(ax2,x,YGiner, 'c', 'linewidth',2);
        hold on
        P6=plot(ax2,x,YPuck, '-.m', 'linewidth',2);
        hold on
        if strcmp(podaci, 'simulirano')
            P5=plot(ax2,xpodaci, zaprikazOriginala, ':', 'color',[0.5 0.5 0.5],'linewidth',2);
            legend(ax2,[ P2 P3 P5 P6 P1], 'HT_{\it \xi=2}','HT_{Hewitt}','Target','Puck','Data', 'location','best','FontSize',FS)
        elseif strcmp(podaci, 'uvezeno')
            P5=plot(ax2,Vv_provjera, zaprikazOriginala_RJ, ':', 'color',[0.5 0.5 0.5],'linewidth',2);
            legend(ax2,[ P2 P3 P5 P6 P1], 'HT_{\it \xi=2}','HT_{Hewitt}','Target','Puck','Data', 'location','best','FontSize',FS)
        else
            P5=plot(ax2,njihovFitVv,njihovFit, ':', 'color',[0.5 0.5 0.5],'linewidth',2);
            legend(ax2,[ P2 P3 P5 P6 P1], 'HT_{\it \xi=2}','HT_{Hewitt}','HT_{[4]}','Puck','Data', 'location','best','FontSize',FS)
        end
        legend boxoff
        xlabel(['Fiber volume','\it V_f '],'FontSize',20);ylabel(['Elastic modulus','{\it E_2}', ' [GPa]'],'FontSize',20)
        grid on
        
        E2v_kv=fitresultKvadrat.E2v;
        Em_kv=fitresultKvadrat.Em;
        E2v_hw=fitresultHewit.E2v;
        Em_hw=fitresultHewit.Em;
        E2v_pk=fitresultPuck.E2v;
        Em_pk=fitresultPuck.Em;
        %         E2v_gn=fitresultGiner.E2v;
        %         Em_gn=fitresultGiner.Em;
        
        text(ax2,0.36,28,['HT-original: E_{2f}=', num2str(E2v_kv),' E_m=', num2str(Em_kv)])
        text(ax2,0.36,26,['Hewit: E_{2f}=', num2str(E2v_hw),' E_m=', num2str(Em_hw)])
        %         text(ax2,0.2,24,['Giner: E_{2f}=', num2str(E2v_gn),' E_m=', num2str(Em_gn)])
        text(ax2,0.36,24,['Puck: E_{2f}=', num2str(E2v_pk),' E_m=', num2str(Em_pk)])
        if strcmp(podaci, 'simulirano')
            text(ax2, 0.36,22,['Terget : E_{2f}=',num2str(E2v_provjera),' E_m=',num2str(Em_provjera)])
        elseif strcmp(podaci, 'uvezeno')
            text(ax2, 0.36,22,'Terget : E_{2f}=65 E_m=4')
        else
            text(ax2, 0.36,22,'Originalno_{clanak} : E_{2f}=73.1 E_m=3.45')
        end
        
        % ci = confint(fitresult,0.95)
        xlim(ax2,[0.35,0.81])
end

%  daci=[xpodaci' ypodaci'];
%  filename = 'podaciZaRJ_5.xlsx';
% %
%  xlswrite(filename,daci,'Sheet1','A1')
% napomena:  podaci do  10 su iz MCMC za RJ
% podaciZaRJ_1 model original mali šumiæ   0.1     E2v_provjera=70; Em_provjera=4; podataka 12   seed=100;
% podaciZaRJ_2 model original  šumiæ   1     E2v_provjera=70; Em_provjera=4; podataka 15   seed=100;
% podaciZaRJ_3 model original  šumiæ   1.2     E2v_provjera=70; Em_provjera=4; podataka 12   seed=3;
% podaciZaRJ_4 model original  šumiæ   1     E2v_provjera=70; Em_provjera=4; podataka 12   seed=3;
% podaciZaRJ_5 model original  šumiæ   1     E2v_provjera=65; Em_provjera=4; podataka 12   seed=3;


