clear all ;clc ;close all
podaciExp =importdata('E2podaci.dat');
kojePodatkeKoristis='experiment' % experiment simulirano uvezeno
sum=8; % kada je šum mali onda jumpPr može biti manji i dobro æe raditi!
dumiModel='Hewitt';%original Hewitt  Puck;
jumpPr = 0.9; % probability of jump
brPodataka=12;

switch kojePodatkeKoristis
    case 'experiment'
        xpodaci=podaciExp(:,1)';
        ypodaci=podaciExp(:,2)';
    case 'simulirano'
        E2v_provjera=65;
        Em_provjera=4;
        Vv_provjera=linspace(0.3,0.8,brPodataka);
        if strcmp(dumiModel,'Puck') % kvadratni raspored kruznih vlakana
            E2_provjera=E_2Puck(Vv_provjera,E2v_provjera,Em_provjera); % Puck
        else
            E2_provjera=E_2HalpinTsai(Vv_provjera,E2v_provjera,Em_provjera,dumiModel); %kvadrat Hewitt Giner
        end
        seed=1; rand( 'state' , seed ); randn('state',seed );
        E2_provjera=E2_provjera+normrnd(0,sum,size(E2_provjera));
        xpodaci=Vv_provjera;
        ypodaci=E2_provjera;
    case 'uvezeno'
        % 4 odgovara MCMC sa više šuma 3 odgovara MCMC sa manje šuma
        % ogovara MCMC sa prve slike!
        nodes= xlsread('podaciZaRJ_10.xlsx');
        xpodaci=nodes(:,1)';
        ypodaci=nodes(:,2)';
 
end

% plot(xpodaci,ypodaci)

%%
FS=18;
Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['data' ],...
    'NumberTitle','off','Position',[0.1 0.1 0.8 0.72]);
ax = axes('Parent',Fig_Optimist,'Position',[0.1 0.15 0.83 0.75],...
    'FontName','Times New Roman','FontSize',FS,...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5);hold on,
plot(ax,xpodaci,ypodaci, 'sr')
xlabel(' Fiber volume V_v','FontSize',FS);ylabel(' Elastic modulus E_2  [GPa]','FontSize',FS)
grid on

%% za raèun prikaz
E2v_linspace=linspace(30,100,1e3);
Em_linspace=linspace(1,10,1e3);
s_raspon=linspace(1e-3,5,1e3);
%% parametri za Normal proposal distribuciju


switch kojePodatkeKoristis
    case 'experiment'
        proposal_E2v=1;
        proposal_Em=0.1;
        proposal_s=0.1;
        %%priori
        a_1= 73; b_1=5 ;
        a_2= 3; b_2=1 ;
        a_3= 1.5; b_3=1 ;
    case 'simulirano'
        %         proposal_E2v=0.5;
        %         proposal_Em=0.1;
        %         proposal_s=0.01;
        %         %%priori
        %         a_1= 60; b_1=5 ;
        %         a_2= 3; b_2=0.5 ;
        %         a_3= 0.5; b_3=0.5 ;
        % za puno puno šuma

%         %% za sliku 7
%         proposal_E2v=1;
%         proposal_Em=0.1;
%         proposal_s=0.1;
%         a_1= 55; b_1=10 ;
%         a_2= 3; b_2=1 ;
%         a_3= 0.5; b_3=0.1 ;
% %%priori za sliku 8
        proposal_E2v=4;
        proposal_Em=1;
        proposal_s=0.5;
        a_1= 55; b_1=10 ;
        a_2= 3; b_2=1 ;
        a_3= 3; b_3=0.5 ;
    case 'uvezeno'
        %     % za NORMAL!!!!!!
        a_1= 50; b_1=10 ; % kada je mali šum prior oupæe nije bitan!
        a_2= 3; b_2=0.5 ;
        a_3= 1; b_3=0.5 ;
        %parametri za Normal proposal distribuciju
        proposal_E2v=1;
        proposal_Em=0.5;
        proposal_s=0.1;
        
end

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

%%
nSamples =numel( xpodaci);

nMCMCSamples = 50000;

modelSamples = zeros(nMCMCSamples, 1);

E2vSamples = zeros(nMCMCSamples, 1);
EmSamples = zeros(nMCMCSamples, 1);
sSamples = zeros(nMCMCSamples, 1);


sampleJ = 1;

kvadrat=1;
Hewitt=2;
Puck=3;

%
modelSamples(1) = [floor(rand*3)+1]; % random initialization / SLUÈAJNO BIRA JEDNOG OD 3 modela sa kojim æe poèeti
%  modelSamples(1) =1;
seed=50; rand( 'state' , seed ); randn('state',seed ); % see je bija 50 bitno za eksperimantalne podatke 
E2v_0 = unifrnd( min(E2v_linspace) , max(E2v_linspace));
Em_0= unifrnd( min(Em_linspace) , max(Em_linspace));
s_0 = unifrnd( min(s_raspon) , max(s_raspon));

E2vSamples(1) = E2v_0;
EmSamples(1) = Em_0;
sSamples(1) = s_0;

%% simulation parameters


prop=[proposal_E2v,proposal_Em,proposal_s];

nJump = 0;
nUpdate = 0;
nTryJump = 0;
nTryUpdate = 0;

seed=20; rand( 'state' , seed ); randn('state',seed );% bija je 20, bitno za ekserimentalne podatke!
% priori modela
priorHewitt=1;
priorKvadrat=1;
priorPuck=1;
tic
while 1
    
    if sampleJ == nMCMCSamples, break; end
    
    u = rand;
    
    model = modelSamples(sampleJ);
    
    E2v0=E2vSamples(sampleJ);
    Em0= EmSamples(sampleJ);
    s0=sSamples(sampleJ);
    
    
    if u<jumpPr % jump
        
        switch( model )
            case kvadrat
                
                pmodel = [Hewitt  Puck];
                
                % kvadrat=1;
                % Hewitt=2;
                % Puck=3;
                
                pc0 = E2v0;
                pc1 = Em0;
                pc2 = s0;
                J = [ 1 0 0 ; 0 1 0 ;0 0 1];
                
                A1=exp(likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'Hewitt')-likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'original'))*priorHewitt/priorKvadrat*(1-jumpPr)/jumpPr*det(J);
                A2=exp(likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'Puck')-likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'original'))*priorPuck/priorKvadrat*(1-jumpPr)/jumpPr*det(J);
                
                [Amax, locAmax ]=max([A1,A2]);
                A = min([1, Amax]);
                
            case Hewitt
                
                pmodel = [kvadrat  Puck ];
                
                pc0 = E2v0;
                pc1 = Em0;
                pc2 = s0;
                J = [ 1 0 0; 0 1 0;0 0 1 ];
                
                A1=exp(likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'original')-likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'Hewitt'))*priorKvadrat/priorHewitt*(1-jumpPr)/jumpPr*det(J);
                A2=exp(likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'Puck')-likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'Hewitt'))*priorPuck/priorHewitt*(1-jumpPr)/jumpPr*det(J);
                
                [Amax, locAmax ]=max([A1,A2]);
                A = min([1, Amax]);
                
                
            case Puck
                
                pmodel = [kvadrat Hewitt];
                
                pc0 = E2v0;
                pc1 = Em0;
                pc2 = s0;
                J = [ 1 0 0; 0 1 0;0 0 1 ];
                
                A1=exp(likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'original')-likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'Puck'))*priorKvadrat/priorPuck*(1-jumpPr)/jumpPr*det(J);
                A2=exp(likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'Hewitt')-likelihoodHTrjmcmc(xpodaci,ypodaci,pc0,pc1,pc2,'Puck'))*priorHewitt/priorPuck*(1-jumpPr)/jumpPr*det(J);
                
                [Amax, locAmax ]=max([A1,A2]);
                A = min([1, Amax]);
                
            otherwise
                sprintf('unknown model %i', model);
        end
        
        if rand<A % accept jump
            nJump = nJump + 1;
            sampleJ = sampleJ + 1;
            %             modelSamples(sampleJ) = randsample(pmodel,1);% bilo ovako modelSamples(sampleJ) =pmodel;
            modelSamples(sampleJ) = pmodel(locAmax);
            E2vSamples(sampleJ) = pc0;
            EmSamples(sampleJ) = pc1;
            sSamples(sampleJ) = pc2;
        end
        nTryJump = nTryJump + 1;
        
        %fprintf('now %i, A=%f\n', model, A);
        
    else % update
        
        switch( model )
            case kvadrat
                %  sampleHTrjmcmc(x, y, poc,prop,pri, modelKompozita,nSamples)
                samp = sampleHTrjmcmc(xpodaci,ypodaci, [E2v0,Em0,s0],prop,pri, 'original');
            case Hewitt
                samp = sampleHTrjmcmc(xpodaci,ypodaci, [E2v0,Em0,s0],prop, pri,'Hewitt');
                %             case Giner
                %                 samp = sampleHTrjmcmc(xpodaci,ypodaci, [E2v0,Em0,s0],prop,pri, 'Giner');
            case Puck
                samp = sampleHTrjmcmc(xpodaci,ypodaci, [E2v0,Em0,s0],prop,pri, 'Puck');
            otherwise
                sprintf('unknown model %i', model);
        end
        
        nUpdate = nUpdate + 1;
        sampleJ = sampleJ + 1;
        modelSamples(sampleJ) = model;
        E2vSamples(sampleJ) = samp(1);
        EmSamples(sampleJ) =samp(2);
        sSamples(sampleJ) = samp(3);
        
    end
    
end

%%

counts = histc(modelSamples, 1:3 );
%%
vjerojatnostiModela=[round(counts(1)/nMCMCSamples*100,1),round(counts(2)/nMCMCSamples*100,1),round(counts(3)/nMCMCSamples*100,1)];
toc
% Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['sve' ],...
%     'NumberTitle','off','Position',[0.01 0.1 0.45 0.72]);
% ax = axes('Parent',Fig_Optimist,'Position',[0.1 0.15 0.8 0.75],...
%     'FontName','Times New Roman','FontSize',12,...
%     'TickDir','in','TickLength',[0.007 0.007],'XTickLabel',...
%     {'','',['M_1(',num2str(vjerojatnostiModela(1)),')%'],...
%     '', ['  M_2(',num2str(vjerojatnostiModela(2)),')%'],...
%     '', ['    M_3(',num2str(vjerojatnostiModela(3)),')%']},'XTickLabelRotation',0,'LineWidth',1.5);hold on,
% barHMCMC = bar(ax,counts/nMCMCSamples);
% txt1=[ 'M_1 original; M_2 Hewitt ;  M_3 Puck'];
% text(0.0, max(counts/nMCMCSamples)*1.05, txt1, 'fontname','times new roman', 'fontsize', 12)
%%

kvadratSampleIed = find(modelSamples==kvadrat);
HewittSampleIed = find(modelSamples==Hewitt);
% GinerSampleIed = find(modelSamples==Giner);
PuckSampleIed = find(modelSamples==Puck);

aa=[0:numel(kvadratSampleIed)-1];
bb=[0:numel(HewittSampleIed)-1];
% cc=[0:numel(GinerSampleIed)-1];
dd=[0:numel(PuckSampleIed)-1];

bIn_1=floor(numel(kvadratSampleIed)*0.1);
bIn_2=floor(numel(HewittSampleIed)*0.1);
% bIn_3=floor(numel(GinerSampleIed)*0.1);
bIn_3=floor(numel(PuckSampleIed)*0.1);

% srednjVrijednost=[ mean(E2vSamples(kvadratSampleIed(bIn_1:end))) mean(EmSamples(kvadratSampleIed(bIn_1:end))) mean(sSamples(kvadratSampleIed(bIn_1:end)));
%     mean(E2vSamples(HewittSampleIed(bIn_2:end))) mean(EmSamples(HewittSampleIed(bIn_2:end))) mean(sSamples(HewittSampleIed(bIn_2:end))) ;
%     mean(E2vSamples(GinerSampleIed(bIn_3:end))) mean(EmSamples(GinerSampleIed(bIn_3:end))) mean(sSamples(GinerSampleIed(bIn_3:end))) ;
%     mean(E2vSamples(PuckSampleIed(bIn_3:end))) mean(EmSamples(PuckSampleIed(bIn_3:end))) mean(sSamples(PuckSampleIed(bIn_3:end)))];

srednjVrijednost=[ mean(E2vSamples(kvadratSampleIed(bIn_1:end))) mean(EmSamples(kvadratSampleIed(bIn_1:end))) mean(sSamples(kvadratSampleIed(bIn_1:end)));
    mean(E2vSamples(HewittSampleIed(bIn_2:end))) mean(EmSamples(HewittSampleIed(bIn_2:end))) mean(sSamples(HewittSampleIed(bIn_2:end))) ;
    mean(E2vSamples(PuckSampleIed(bIn_3:end))) mean(EmSamples(PuckSampleIed(bIn_3:end))) mean(sSamples(PuckSampleIed(bIn_3:end)))];

%%
Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['MCMC' ],...
    'NumberTitle','off','Position',[0.1 0.05 0.8 0.85]);
ax1 = axes('Parent',Fig_Optimist,'Position',[0.1 0.75 0.8 0.23],...
    'FontName','Times New Roman','FontSize',14,...
    'TickDir','in','TickLength',[0.007 0.007]);hold on;
ax2 = axes('Parent',Fig_Optimist,'Position',[0.1 0.4 0.8 0.23],...
    'FontName','Times New Roman','FontSize',14,...
    'TickDir','in','TickLength',[0.007 0.007]);hold on;
ax3 = axes('Parent',Fig_Optimist,'Position',[0.1 0.1 0.8 0.23],...
    'FontName','Times New Roman','FontSize',14,...
    'TickDir','in','TickLength',[0.007 0.007]);hold on;
% ax4 = axes('Parent',Fig_Optimist,'Position',[0.1 0.15 0.8 0.15],...
%     'FontName','Times New Roman','FontSize',14,...
%     'TickDir','in','TickLength',[0.007 0.007]);hold on;

plot(ax1,aa,E2vSamples(kvadratSampleIed),'xb', 'markersize',2,'linewidth',0.5);
hold on
plot(ax1,aa,EmSamples(kvadratSampleIed),'xm', 'markersize',2,'linewidth',0.5);
hold on
plot(ax1,aa,sSamples(kvadratSampleIed),'xg', 'markersize',2,'linewidth',0.5);

% srednjVrijednost(1,1) srednjVrijednost(1,2) srednjVrijednost(1,3)
% srednjVrijednost(2,1) srednjVrijednost(2,2) srednjVrijednost(2,3)
% srednjVrijednost(3,1)srednjVrijednost(3,2) srednjVrijednost(3,3)
hold on
plot(ax1,[0 max(aa)],[srednjVrijednost(1,1) srednjVrijednost(1,1) ],'-b', 'linewidth',0.5);
hold on
plot(ax1,[0 max(aa)],[srednjVrijednost(1,2) srednjVrijednost(1,2)],'-m', 'linewidth',0.5);
hold on
plot(ax1,[0 max(aa)],[srednjVrijednost(1,3) srednjVrijednost(1,3)],'-g', 'linewidth',0.5);
ylabel(ax1,' original','FontSize',14)

plot(ax2,bb,E2vSamples(HewittSampleIed),'xb', 'markersize',2,'linewidth',0.5);
hold on
plot(ax2,bb,EmSamples(HewittSampleIed),'xm', 'markersize',2,'linewidth',0.5);
hold on
plot(ax2,bb,sSamples(HewittSampleIed),'xg', 'markersize',2,'linewidth',0.5);
hold on
plot(ax2,[0 max(bb)],[srednjVrijednost(2,1) srednjVrijednost(2,1) ],'-b', 'linewidth',0.5);
hold on
plot(ax2,[0 max(bb)],[ srednjVrijednost(2,2)  srednjVrijednost(2,2)],'-m', 'linewidth',0.5);
hold on
plot(ax2,[0 max(bb)],[srednjVrijednost(2,3) srednjVrijednost(2,3)],'-g', 'linewidth',0.5)
ylabel(ax2,' Hewitt','FontSize',14)

% plot(ax4,cc,E2vSamples(GinerSampleIed),'xb', 'markersize',2,'linewidth',0.5);
% hold on
% plot(ax4,cc,EmSamples(GinerSampleIed),'xm', 'markersize',2,'linewidth',0.5);
% hold on
% plot(ax4,cc,sSamples(GinerSampleIed),'xg', 'markersize',2,'linewidth',0.5);
% hold on
% plot(ax4,[0 max(cc)],[srednjVrijednost(3,1) srednjVrijednost(3,1)],'-b', 'linewidth',0.5);
% hold on
% plot(ax4,[0 max(cc)],[srednjVrijednost(3,2) srednjVrijednost(3,2)],'-m', 'linewidth',0.5);
% hold on
% plot(ax4,[0 max(cc)],[srednjVrijednost(3,3) srednjVrijednost(3,3)],'-g', 'linewidth',0.5);
% xlabel(ax4,' uzorci','FontSize',14); ylabel(ax3,' Giner','FontSize',14)

plot(ax3,dd,E2vSamples(PuckSampleIed),'xb', 'markersize',2,'linewidth',0.5);
hold on
plot(ax3,dd,EmSamples(PuckSampleIed),'xm', 'markersize',2,'linewidth',0.5);
hold on
plot(ax3,dd,sSamples(PuckSampleIed),'xg', 'markersize',2,'linewidth',0.5);
hold on
plot(ax3,[0 max(dd)],[srednjVrijednost(3,1) srednjVrijednost(3,1)],'-b', 'linewidth',0.5);
hold on
plot(ax3,[0 max(dd)],[srednjVrijednost(3,2) srednjVrijednost(3,2)],'-m', 'linewidth',0.5);
hold on
plot(ax3,[0 max(dd)],[srednjVrijednost(3,3) srednjVrijednost(3,3)],'-g', 'linewidth',0.5);
xlabel(ax3,' uzorci','FontSize',14); ylabel(ax3,' Puck','FontSize',14)

text(ax1,max(aa)*0.05,20,['E_{2v}=',num2str(srednjVrijednost(1,1))], 'Color','b','FontSize',12);
text(ax1,max(aa)*0.2,20,['E_{m}=',num2str(srednjVrijednost(1,2))], 'Color','m','FontSize',12);
text(ax1,max(aa)*0.35,20,['\sigma=',num2str(srednjVrijednost(1,3))], 'Color','g','FontSize',12);

text(ax2,max(bb)*0.05,20,['E_{2v}=',num2str(srednjVrijednost(2,1))], 'Color','b','FontSize',12);
text(ax2,max(bb)*0.2,20,['E_{m}=',num2str(srednjVrijednost(2,2))], 'Color','m','FontSize',12);
text(ax2,max(bb)*0.35,20,['\sigma=',num2str(srednjVrijednost(2,3))], 'Color','g','FontSize',12);
%
% text(ax3,max(cc)*0.05,20,['E_{2v}=',num2str(srednjVrijednost(3,1))], 'Color','b','FontSize',12);
% text(ax3,max(cc)*0.2,20,['E_{m}=',num2str(srednjVrijednost(3,2))], 'Color','m','FontSize',12);
% text(ax3,max(cc)*0.35,20,['\sigma=',num2str(srednjVrijednost(3,3))], 'Color','g','FontSize',12);

text(ax3,max(dd)*0.05,20,['E_{2v}=',num2str(srednjVrijednost(3,1))], 'Color','b','FontSize',12);
text(ax3,max(dd)*0.2,20,['E_{m}=',num2str(srednjVrijednost(3,2))], 'Color','m','FontSize',12);
text(ax3,max(dd)*0.35,20,['\sigma=',num2str(srednjVrijednost(3,3))], 'Color','g','FontSize',12);

%%

Vv_model=0.2:0.01:0.8;
E2v_Bayes_M1=srednjVrijednost(1,1);
Em_Bayes_M1=srednjVrijednost(1,2);
E2_model_M1=E_2HalpinTsai(Vv_model,E2v_Bayes_M1,Em_Bayes_M1,'original');

E2v_Bayes_M2=srednjVrijednost(2,1);
Em_Bayes_M2=srednjVrijednost(2,2);
E2_model_M2=E_2HalpinTsai(Vv_model,E2v_Bayes_M2,Em_Bayes_M2,'Hewitt');

E2v_Bayes_M3=srednjVrijednost(3,1);
Em_Bayes_M3=srednjVrijednost(3,2);
E2_model_M3=E_2Puck(Vv_model,E2v_Bayes_M3,Em_Bayes_M3);

% E2v_Bayes_M4=srednjVrijednost(4,1)
% Em_Bayes_M4=srednjVrijednost(4,2)
% E2_model_M4=E_2Puck(Vv_model,E2v_Bayes_M4,Em_Bayes_M4);


%{
FS=18;
Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['fitano' ],...
    'NumberTitle','off','Position',[0.47 0.1 0.51 0.72]);
ax = axes('Parent',Fig_Optimist,'Position',[0.13 0.15 0.83 0.75],...
    'FontName','Times New Roman','FontSize',FS,...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5);hold on,
[~, w]=max(vjerojatnostiModela);
if w==1
    P1=plot(ax,Vv_model,E2_model_M1, '-k', 'linewidth',3);
    hold on
    P2=plot(ax,Vv_model,E2_model_M2, '--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    P3=plot(ax,Vv_model,E2_model_M3, ':','color',[0.5 0.5 0.5] ,'linewidth',1);
    %     hold on
    %     P4=plot(ax,Vv_model,E2_model_M4, '-.','color',[0.5 0.5 0.5], 'linewidth',1);
    hold on
    P5=plot(ax,xpodaci,ypodaci, 'sr');
elseif w==2
    P1=plot(ax,Vv_model,E2_model_M1, '--' ,'color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    P2=plot(ax,Vv_model,E2_model_M2, '-b', 'linewidth',3);
    hold on
    P3=plot(ax,Vv_model,E2_model_M3, ':','color',[0.5 0.5 0.5], 'linewidth',1);
    %     hold on
    %     P4=plot(ax,Vv_model,E2_model_M4, '-.','color',[0.5 0.5 0.5], 'linewidth',1);
    hold on
    P5=plot(ax,xpodaci,ypodaci, 'sr');
elseif w==3
    P1=plot(ax,Vv_model,E2_model_M1, '--','color',[0.5 0.5 0.5], 'linewidth',1);
    hold on
    P2=plot(ax,Vv_model,E2_model_M2, ':','color',[0.5 0.5 0.5], 'linewidth',1);
    hold on
    P3=plot(ax,Vv_model,E2_model_M3, 'c', 'linewidth',3);
    %     hold on
    %     P4=plot(ax,Vv_model,E2_model_M4, '-.','color',[0.5 0.5 0.5], 'linewidth',1);
    hold on
    P5=plot(ax,xpodaci,ypodaci, 'sr');
    % elseif w==4
    %     P1=plot(ax,Vv_model,E2_model_M1, '--','color',[0.5 0.5 0.5], 'linewidth',1);
    %     hold on
    %     P2=plot(ax,Vv_model,E2_model_M2, ':b','color',[0.5 0.5 0.5], 'linewidth',1);
    %     hold on
    %     P3=plot(ax,Vv_model,E2_model_M3, '-.','color',[0.5 0.5 0.5], 'linewidth',1);
    %     hold on
    %     P4=plot(ax,Vv_model,E2_model_M4, 'm', 'linewidth',3);
    %     hold on,
    %     P5=plot(ax,xpodaci,ypodaci, 'sr');
end
legend(ax, [P1 P2 P3 P5], ['original ',num2str(vjerojatnostiModela(1)),'%'],['Hewitt ',num2str(vjerojatnostiModela(2)),'%'],...
    ['Puck ',num2str(vjerojatnostiModela(3)),'%'],'Data','location','best',  'FontSize',12)
legend boxoff
xlabel(['Fiber volume','\it V_f '],'FontSize',FS);ylabel(['Elastic modulus','{\it E_2}', ' [GPa]'],'FontSize',FS)
grid on
%}
%% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

Fig_zadnja = figure('units','normalized','Color',[1 1 1],'Name',['izbor modela' ],...
    'NumberTitle','off','Position',[0.1 0.05 0.8 0.85]);
ax01 = axes('Parent',Fig_zadnja,'Position',[0.055 0.15 0.35 0.75],...
    'FontName','Times New Roman','FontSize',20,...
    'TickDir','in','TickLength',[0.007 0.007],'XTickLabel',...
    {' ','M_1', '   M_2', '    M_3'},'XTickLabelRotation',0,'LineWidth',1.5);hold on,
ax02 = axes('Parent',Fig_zadnja,'Position',[0.45 0.15 0.52 0.75],...
    'FontName','Times New Roman','FontSize',20,...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5);hold on,
barHMCMC = bar(ax01,counts/nMCMCSamples);
txt1=[ 'M_1, HT_{\xi=2}    ',num2str(vjerojatnostiModela(1)),'%'];
txt2=['M_2,  HT_{Hewitt}',num2str(vjerojatnostiModela(2)),'%'];
txt3=['M_3,  Puck     ',num2str(vjerojatnostiModela(3)),'%'];
text(ax01,1.7, max(counts/nMCMCSamples), txt1, 'fontname','times new roman', 'fontsize', 15)
text(ax01,1.7, max(counts/nMCMCSamples)*0.93, txt2, 'fontname','times new roman', 'fontsize', 15)
text(ax01,1.7, max(counts/nMCMCSamples)*0.86, txt3, 'fontname','times new roman', 'fontsize', 15)
[~, w]=max(vjerojatnostiModela);
if w==1
    P1=plot(ax02,Vv_model,E2_model_M1, '-k', 'linewidth',3);
    hold on
    P2=plot(ax02,Vv_model,E2_model_M2, '--','color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    P3=plot(ax02,Vv_model,E2_model_M3, ':','color',[0.5 0.5 0.5] ,'linewidth',1);
    hold on
    P5=plot(ax02,xpodaci,ypodaci, 'sr');
elseif w==2
    P1=plot(ax02,Vv_model,E2_model_M1, '--' ,'color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    P2=plot(ax02,Vv_model,E2_model_M2, '-b', 'linewidth',3);
    hold on
    P3=plot(ax02,Vv_model,E2_model_M3, ':','color',[0.5 0.5 0.5], 'linewidth',1);
    hold on
    P5=plot(ax02,xpodaci,ypodaci, 'sr');
elseif w==3
    P1=plot(ax02,Vv_model,E2_model_M1, '--','color',[0.5 0.5 0.5], 'linewidth',1);
    hold on
    P2=plot(ax02,Vv_model,E2_model_M2, ':','color',[0.5 0.5 0.5], 'linewidth',1);
    hold on
    P3=plot(ax02,Vv_model,E2_model_M3, 'c', 'linewidth',3);
    hold on
    P5=plot(ax02,xpodaci,ypodaci, 'sr');
end
legend(ax02, [P1 P2 P3 P5], 'HT_{\xi=2}','HT_{Hewitt} ',...
   'Puck ','Data','location','best',  'FontSize',18)
legend boxoff
xlabel(ax02,['Fiber volume','\it V_f '],'FontSize',20);ylabel(ax02,['Elastic modulus','{\it E_2}', ' [GPa]'],'FontSize',20)
grid on
xlim(ax02,[0.3 0.8])
annotation(Fig_zadnja,'textbox',[0.01 0.92 0.03 0.07],...
    'String','a)',...
    'FontSize',28,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
annotation(Fig_zadnja,'textbox',[0.38 0.92 0.03 0.07],...
    'String','b)',...
    'FontSize',28,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
%% 
return
  daci=[xpodaci' ypodaci'];
  filename = 'podaciZaRJ_13.xlsx';
 xlswrite(filename,daci,'Sheet1','A1')
% napomena: podaci od 10 pa nadalje su iz JR za MCMC
% podaciZaRJ_10 model original, šum 1.2     E2v_provjera=65; Em_provjera=4; podataka 12   seed=3; 
% podaciZaRJ_11 model Hewitt, šum 2     E2v_provjera=65; Em_provjera=4; podataka 12   seed=2; 
% podaciZaRJ_12 model Hewitt, šum 10     E2v_provjera=65; Em_provjera=4; podataka 12   seed=2; 
% podaciZaRJ_13 model Hewitt, šum 2     E2v_provjera=65; Em_provjera=4; podataka 12   seed=8;  Vv_provjera=linspace(0.3,0.8,brPodataka);



