function stipePlotHT(N,nbins, samples, estimated,  burnIn, name)
%% za prvi parametar
bins_1_cl= linspace( (min(samples(burnIn:end,1))) , (max(samples(burnIn:end,1))) , nbins*10 );
counts_1_cl = hist( samples(burnIn:end,1) , bins_1_cl );
countsNorm_1_cl=counts_1_cl/sum(counts_1_cl);

% LowIndx_1=find(cumsum(countsNorm_1_cl)<=0.05 );
% UpwIndx_1=find(cumsum(countsNorm_1_cl)>=0.95 );
% binsL_1=bins_1_cl(LowIndx_1(end))
% binsU_1=bins_1_cl(UpwIndx_1(1))

bins_1 = linspace( (min(samples(burnIn:end,1))) , (max(samples(burnIn:end,1))) , nbins );
counts_1 = hist( samples(burnIn:end,1) , bins_1 );
countsNorm_1=counts_1/sum(counts_1);
maksPrvi=max(countsNorm_1/mean(diff(bins_1)));
%%  za drugi parametar
bins_2_cl = linspace( (min(samples(burnIn:end,2))) , (max(samples(burnIn:end,2))) , nbins*10 );
counts_2_cl = hist( samples(burnIn:end,2) , bins_2_cl );
countsNorm_2_cl =counts_2_cl /sum(counts_2_cl );
% 
% LowIndx_2=find(cumsum(countsNorm_2_cl)<=0.05 );
% UpwIndx_2=find(cumsum(countsNorm_2_cl)>=0.95 );
% % binsL_2=bins_2_cl(LowIndx_2(end))
% binsU_2=bins_2_cl(UpwIndx_2(1))

bins_2 = linspace( (min(samples(burnIn:end,2))) , (max(samples(burnIn:end,2))) , nbins );
counts_2 = hist( samples(burnIn:end,2) , bins_2 );
countsNorm_2 =counts_2 /sum(counts_2 );
maksDrugi=max(countsNorm_2 /mean(diff(bins_2)));
%%  za treci parametar
bins_3_cl = linspace( (min(samples(burnIn:end,3))) , (max(samples(burnIn:end,3))) , nbins*10 );
counts_3_cl = hist( samples(burnIn:end,3) , bins_3_cl );
countsNorm_3_cl =counts_3_cl /sum(counts_3_cl );

% LowIndx_3=find(cumsum(countsNorm_3_cl)<=0.05 );
% UpwIndx_3=find(cumsum(countsNorm_3_cl)>=0.95 );
% binsL_3=bins_3_cl(LowIndx_3(end))
% binsU_3=bins_3_cl(UpwIndx_3(1))

bins_3 = linspace( (min(samples(burnIn:end,3))) , (max(samples(:,3))) , nbins );
counts_3 = hist( samples(burnIn:end,3) , bins_3 );
countsNorm_3 =counts_3 /sum(counts_3 );
maksTreci=max(countsNorm_3 /mean(diff(bins_3)));
%% definiranje slike i osiju
xMin_1= min(samples(:,1))*0.9; xMax_1=max(samples(:,1))*1.05;
xMin_2= min(samples(:,2))*0.9; xMax_2=max(samples(:,2))*1.05;
xMin_3= min(samples(:,3))*0.9; xMax_3=max(samples(:,3))*1.05;

Fig_Optimist = figure('units','normalized','Color',[1 1 1],'Name',['MCMC za ', name ],...
    'NumberTitle','off','Position',[0.1 0.1 0.8 0.72]); % ili jednostavno [0 0 1 1] i sve štima ili [0.003 0.061 0.994 0.83]
ax1 = axes('Parent',Fig_Optimist,'Position',[0.075 0.71 0.25 0.25],...
    'FontName','Times New Roman','FontSize',18,...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5,'TickLength',[0.007 0.007],'XTickLabel',{'','','','',''},...
    'Xlim',[xMin_1 , xMax_1],'Ylim',[0, maksPrvi*1.1]);hold on,
ylabel(' Frequency','FontSize',18)

ax2 = axes('Parent',Fig_Optimist,'Position',[0.075 0.15  0.25 0.5],...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5,...
    'FontName','Times New Roman','FontSize',18, ...
    'Xlim',[xMin_1, xMax_1],'Ylim',[-2000, N]);hold on,
ylabel(' Samples','FontSize',18);hold on
xlabel('\it  E_{2f }','FontSize',18)
 
ax3 = axes('Parent',Fig_Optimist,'Position',[0.395 0.71 0.25 0.25],...
    'FontName','Times New Roman','FontSize',18,...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5,'TickLength',[0.007 0.007],'XTickLabel',{'','','','',''}, ...
    'Xlim',[xMin_2, xMax_2],'Ylim',[0, maksDrugi]);hold on,

ax4 = axes('Parent',Fig_Optimist,'Position',[0.395 0.15  0.25 0.5],...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5,...
    'FontName','Times New Roman','FontSize',18, ...
    'Xlim',[xMin_2, xMax_2],'Ylim',[-2000, N]);hold on,
xlabel('\it  E_m','FontSize',18)
% 
 ax5 = axes('Parent',Fig_Optimist,'Position',[0.705 0.71 0.25 0.25],...
    'FontName','Times New Roman','FontSize',18,...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5,'TickLength',[0.007 0.007],'XTickLabel',{'','','','',''},...
    'Xlim',[xMin_3, xMax_3],'Ylim',[0, maksTreci]);hold on,

ax6 = axes('Parent',Fig_Optimist,'Position',[0.705 0.15  0.25 0.5],...
    'TickDir','in','TickLength',[0.007 0.007],'LineWidth',1.5,...
    'FontName','Times New Roman','FontSize',18, ...
    'Xlim',[xMin_3, xMax_3],'Ylim',[-2000, N]);hold on,
xlabel('\it  \sigma_{LH}','FontSize',18);


%% plotanje
%%
p0=bar(ax1, bins_1 , countsNorm_1/mean(diff(bins_1)) ,'FaceColor','b','EdgeColor',[0 0 0.5],'LineWidth',1  );
hold on
p1=plot(ax1,[estimated(1) estimated(1)],[0 maksPrvi], '-','color', 'r', 'linewidth',2);
  hold on
legend(ax1, [p0 p1],'Histogram','E[{\itE_{2f}}] ','location', 'best' )
legend(ax1,'boxoff')
grid(ax1,'on');

 plot(ax2, samples(:,1) , 1:N ,'-o','color','r','markersize',1,'markeredgecolor', 'b');
hold on
plot(ax2, [xMin_1 , xMax_1] , [burnIn burnIn] , '-','color',  [0.3 0.3 0.3],'linewidth',3);
hold on
% plot(ax2, samples(1,1) , 1 ,'o','linewidth',2,'markersize',6,'markerfacecolor','y','markeredgecolor', 'r');
% legend(ax2, 'Burn in', 'Samples','Initial sample', 'location', 'best' )
legend(ax2,'Samples', 'Burn in',  'location', 'best' )
legend(ax2,'boxoff')
grid(ax2,'on');
%%
p0=bar(ax3, bins_2 , countsNorm_2 /mean(diff(bins_2))  ,'FaceColor','b','EdgeColor',[0 0 0.5],'LineWidth',1);
hold on
p1=plot(ax3,[estimated(2) estimated(2)],[0 maksDrugi],  '-','color', 'r', 'linewidth',2);
hold on
legend(ax3, [p0 p1],'Histogram','E[{\itE_m}] ','location', 'best' )
legend(ax3,'boxoff')
grid(ax3,'on');

plot(ax4, samples(:,2) , 1:N ,'-o','color','r','markersize',1,'markeredgecolor', 'b');
hold on
plot(ax4, [xMin_2, xMax_2] , [burnIn burnIn] , '-','color',  [0.3 0.3 0.3],'linewidth',3);

% hold on
% plot(ax4, samples(1,2) , 1 ,'o','linewidth',2,'markersize',6,'markerfacecolor','y','markeredgecolor', 'r');
% legend(ax4, 'Burn in', 'Samples','Initial sample', 'location', 'best' )
legend(ax4, 'Samples','Burn in',  'location', 'best' )
legend(ax4,'boxoff')
grid(ax4,'on');
%%
p0=bar(ax5,bins_3 , countsNorm_3 /mean(diff(bins_3))  ,'FaceColor','b','EdgeColor',[0 0 0.5],'LineWidth',1);
hold on
p1=plot(ax5,[estimated(3) estimated(3)],[0 maksTreci],  '-','color', 'r', 'linewidth',2);
hold on
legend(ax5, [p0 p1],'Histogram','E[\sigma_{LH}] ','location', 'best' )
legend(ax5,'boxoff')
grid(ax5,'on');

plot(ax6, samples(:,3) , 1:N ,'-o','color','r','markersize',1,'markeredgecolor', 'b');
hold on
plot(ax6,  [xMin_3, xMax_3] , [burnIn burnIn] , '-','color',  [0.3 0.3 0.3],'linewidth',3);
% hold on
% plot(ax6, samples(1,3) , 1 ,'o','linewidth',2,'markersize',6,'markerfacecolor','y','markeredgecolor', 'r');
% legend(ax6, 'Burn in', 'Samples','Initial sample', 'location', 'best' )
legend(ax6, 'Samples','Burn in',  'location', 'best' )
legend(ax6,'boxoff')
grid(ax6,'on');
end




