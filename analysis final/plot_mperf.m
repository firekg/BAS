function plot_mperf(PERF, mPERF)
%% Plot: PERF, mPERF

figure(1);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

%plot PERF
pd=0.04; %distance between plots
pws=0.2; %plot width
phs=0.16; %plot height
lm=0.00; %left margin
ls1=lm+2*pd; %1st column left start position
ts2=0.85; %1st row top start position
x=5:5:25;
xm = 1:25;
% PERF has 4 layers: revealing number x5(1:5:25); parti x4; condiitons x3 (AL,PLR,PLB); moments x2 (mean,SD)
for ii=1:4
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts2;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
    plot([0,30],[0.5,0.5],':','Color',[0.5,0.5,0.5]);  hold on;  plot([0,30],[1,1],':','Color',[0.5,0.5,0.5]);
    % mPERF
    plot(xm, mPERF(:,ii,3,1), 'k-');
    a_ErrorShade(xm, mPERF(:,ii,3,1)-mPERF(:,ii,3,2), mPERF(:,ii,3,1)+mPERF(:,ii,3,2), [0,0,0], 1);
    plot(xm, mPERF(:,ii,4,1), '-', 'Color', [0.5, 0.5, 0.5]);
    a_ErrorShade(xm, mPERF(:,ii,4,1)-mPERF(:,ii,4,2), mPERF(:,ii,4,1)+mPERF(:,ii,4,2), 0.5*[1,1,1], 1);
    plot(xm, mPERF(:,ii,2,1), 'b-');
    a_ErrorShade(xm, mPERF(:,ii,2,1)-mPERF(:,ii,2,2), mPERF(:,ii,2,1)+mPERF(:,ii,2,2), [0,0,1], 1);
    plot(xm, mPERF(:,ii,1,1), 'r-');
    a_ErrorShade(xm, mPERF(:,ii,1,1)-mPERF(:,ii,1,2), mPERF(:,ii,1,1)+mPERF(:,ii,1,2), [1,0,0], 1);
    plot(xm, mPERF(:,ii,5,1), 'g-');
    a_ErrorShade(xm, mPERF(:,ii,5,1)-mPERF(:,ii,5,2), mPERF(:,ii,5,1)+mPERF(:,ii,5,2), [0,1,0], 1);    
    % PERF
    errorbar(x,PERF(:,ii,3,1),PERF(:,ii,3,2),'kd','MarkerFaceColor',[0,0,0],'MarkerSize',3);
    errorbar(x,PERF(:,ii,4,1),PERF(:,ii,4,2),'^', 'Color',0.5*[1,1,1],'MarkerFaceColor',0.5*[1,1,1],'MarkerSize',3);
    errorbar(x,PERF(:,ii,2,1),PERF(:,ii,2,2),'bs','MarkerFaceColor',[0,0,1],'MarkerSize',3);
    errorbar(x,PERF(:,ii,1,1),PERF(:,ii,1,2),'ro','MarkerFaceColor',[1,0,0],'MarkerSize',3);    
    ylabel('P(correct)');  % ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);  set(gca,'FontSize',10);
    set(gca,'Ytick',0.3:0.1:1,'YTickLabel',{'','','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.01,0]);
    axis([0 27.5 0.31 1.05]);  box off;
    title(['P',num2str(ii)]);
    set(gca,'FontSize',10);
end

set(gcf,'Color',[1,1,1]);

% % save figure
% export_fig('PerfInfo.pdf')

end