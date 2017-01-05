function plot_perf_info(PERF, INFO, mPERF)
%% Plot: model PERF, INFO, INFO heurisitcs

figure(1);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

% initialize some things
% pd=0.005; %distance between plots
% lm=0.05; %left margin
% ps=0.07; %plot size
% ls1=lm+2*pd; %1st column left start position
% ts1=0.85; %1st row top start position

%plot PERF
pd=0.04; %distance between plots
pws=0.2; %plot width
phs=0.16; %plot height
lm=0.00; %left margin
ls1=lm+2*pd; %1st column left start position
ts2=0.45; %1st row top start position
x=5:5:25;
xm = 1:25;
% PERF has 4 layers: revealing number x5(1:5:25); parti x4; condiitons x3 (AL,PLR,PLB); moments x2 (mean,SD)
for ii=1:4
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts2;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
    plot([0,30],[0.5,0.5],':','Color',[0.5,0.5,0.5]);  hold on;  plot([0,30],[1,1],':','Color',[0.5,0.5,0.5]);
    
    % mPERF
    % noiseless BAS
    plot(xm, mPERF(:,ii,3,1), 'k-');
    a_ErrorShade(xm, mPERF(:,ii,3,1)-mPERF(:,ii,3,2), mPERF(:,ii,3,1)+mPERF(:,ii,3,2), [0,0,0], 1);
    
    % anti-BAS
%     plot(xm, mPERF(:,ii,4,1), '-', 'Color', [0.5, 0.5, 0.5]);
%     a_ErrorShade(xm, mPERF(:,ii,4,1)-mPERF(:,ii,4,2), mPERF(:,ii,4,1)+mPERF(:,ii,4,2), 0.5*[1,1,1], 1);
    
    % random
    plot(xm, mPERF(:,ii,2,1), 'b-');
    a_ErrorShade(xm, mPERF(:,ii,2,1)-mPERF(:,ii,2,2), mPERF(:,ii,2,1)+mPERF(:,ii,2,2), [0,0,1], 1);
    
    % active
    plot(xm, mPERF(:,ii,1,1), 'r-');
    a_ErrorShade(xm, mPERF(:,ii,1,1)-mPERF(:,ii,1,2), mPERF(:,ii,1,1)+mPERF(:,ii,1,2), [1,0,0], 1);
    
%     % noisy BAS
    cl=[0.3,0.5,0.3];
    plot(xm, mPERF(:,ii,5,1), '-', 'Color', cl);
    a_ErrorShade(xm, mPERF(:,ii,5,1)-mPERF(:,ii,5,2), mPERF(:,ii,5,1)+mPERF(:,ii,5,2), cl, 1);    
    % prob-match BAS
    cl=[0.3,0.3,0.5];
    plot(xm, mPERF(:,ii,5,1), '-', 'Color', cl);
    a_ErrorShade(xm, mPERF(:,ii,6,1)-mPERF(:,ii,6,2), mPERF(:,ii,6,1)+mPERF(:,ii,6,2), cl, 1);    
    % PERF
    errorbar(x,PERF(:,ii,3,1),PERF(:,ii,3,2),'kd','MarkerFaceColor',[0,0,0],'MarkerSize',3);
%     errorbar(x,PERF(:,ii,4,1),PERF(:,ii,4,2),'^', 'Color',0.5*[1,1,1],'MarkerFaceColor',0.5*[1,1,1],'MarkerSize',3);
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

%plot INFO
pd=0.04; %distance between plots
pws=0.2; %plot width
phs=0.16; %plot height
lm=0.00; %left margin
ls1=lm+2*pd; %1st column left start position
ts2=0.25; %1st row top start position
x=1:25;
% INFO has 4 layers: revealing number x25 (1:25); parti x3; condiitons x5 (AL,PLR,BAS,nBAS,sBAS); moments x2 (mean,SD)
for ii=1:4
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts2;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
    %plot([0,30],[1,1],':','Color',[0.5,0.5,0.5]);   %line on 1
    %BAS->mBAS->PLR->AL->sBAS
    hold on;
%     cl=[0.0,0.0,0.0]; a_ErrorShade(x,INFO(:,ii,3,1)-INFO(:,ii,3,2),INFO(:,ii,3,1)+INFO(:,ii,3,2),cl,1);
%     cl=[0.0,0.0,1.0]; a_ErrorShade(x,INFO(:,ii,2,1)-INFO(:,ii,2,2),INFO(:,ii,2,1)+INFO(:,ii,2,2),cl,1);
%     cl=[1.0,0.0,0.0]; a_ErrorShade(x,INFO(:,ii,1,1)-INFO(:,ii,1,2),INFO(:,ii,1,1)+INFO(:,ii,1,2),cl,1);    
%     cl=[0.3,0.5,0.3]; a_ErrorShade(x,INFO(:,ii,5,1)-INFO(:,ii,5,2),INFO(:,ii,5,1)+INFO(:,ii,5,2),cl,1);
%     cl=[0.3,0.3,0.5]; a_ErrorShade(x,INFO(:,ii,6,1)-INFO(:,ii,6,2),INFO(:,ii,6,1)+INFO(:,ii,6,2),cl,1);
%     cl=[0.0,1.0,0.0]; a_ErrorShade(x,INFO(:,ii,10,1)-INFO(:,ii,10,2),INFO(:,ii,10,1)+INFO(:,ii,10,2),cl,1);

    cl=[0.0,0.0,0.0];  plot(x,INFO(:,ii,3,1),'-','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.0,0.0,1.0];  plot(x,INFO(:,ii,2,1),'-','Color',cl, 'MarkerFaceColor',cl);    
    cl=[1.0,0.0,0.0];  plot(x,INFO(:,ii,1,1),'-','Color',cl, 'MarkerFaceColor',cl);
%     cl=[0.3,0.5,0.3];  plot(x,INFO(:,ii,5,1),'-','Color',cl, 'MarkerFaceColor',cl);
%     cl=[0.3,0.3,0.5];  plot(x,INFO(:,ii,6,1),'-','Color',cl, 'MarkerFaceColor',cl);
%     cl=[0.0,1.0,0.0];  plot(x,INFO(:,ii,10,1),'-','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.5,0.5,0.5];  plot(x,INFO(:,ii,4,1),'-','Color',cl, 'MarkerFaceColor',cl);
    
%     if(ii==1)
%         set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','','','','','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
%         ylabel('Information gain');  % ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);  set(gca,'FontSize',10);
%     else
%         set(gca,'Ytick',0:0.1:1,'YTickLabel',{'','','','','','','','','','',''}, 'TickLength',[0.01,0]);    
%     end
    ylabel('Information (bit)');
    set(gca,'Ytick',0:0.2:1,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'}, 'TickLength',[0.01,0]);    
%     set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
%     set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','','','','','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.01,0]);
%     if (ii==1 | ii==2)
%         axis([0 27.5 0 0.4]);  box off;  
%     elseif (ii==3)
%         axis([0 27.5 0 0.54]);  box off;
%     else
%         axis([0 27.5 0 0.4]);  box off;
%     end
    axis([0 27.5 0 1]);  box off;  % for plotting measured noise
%     axis([0 27.5 0 0.55]);  box off;  % for plotting fitted noise
%     axis([0 27.5 0 0.31]);  box off;  % for plotting heuristics
    set(gca,'FontSize',10);
end
%leg=legend('BAS', 'biased nBAS', 'free-scan', 'random', 'Location','SouthEast'); 
%set(leg,'Color','w');


%plot INFO heuristics
pd=0.04; %distance between plots
pws=0.2; %plot width
phs=0.16; %plot height
lm=0.00; %left margin
ls1=lm+2*pd; %1st column left start position
ts2=0.05; %1st row top start position
x=1:25;
% INFO has 4 layers: revealing number x25 (1:25); parti x3; condiitons x5 (AL,PLR,BAS,nBAS,sBAS); moments x2 (mean,SD)
for ii=1:4
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts2;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
    hold on;
    cl=[1.0,0.0,0.0]; a_ErrorShade(x,INFO(:,ii,1,1)-INFO(:,ii,1,2),INFO(:,ii,1,1)+INFO(:,ii,1,2),cl,1);    
    cl=[1.0,0.4,0.0]; a_ErrorShade(x,INFO(:,ii,7,1)-INFO(:,ii,7,2),INFO(:,ii,7,1)+INFO(:,ii,7,2),cl,1);
    cl=[0.6,0.1,0.9]; a_ErrorShade(x,INFO(:,ii,8,1)-INFO(:,ii,8,2),INFO(:,ii,8,1)+INFO(:,ii,8,2),cl,1);
    cl=[0.6,0.3,0.0]; a_ErrorShade(x,INFO(:,ii,9,1)-INFO(:,ii,9,2),INFO(:,ii,9,1)+INFO(:,ii,9,2),cl,1);

    cl=[1.0,0.0,0.0];  plot(x,INFO(:,ii,1,1),'-','Color',cl, 'MarkerFaceColor',cl);
    cl=[1.0,0.4,0.0];  plot(x,INFO(:,ii,7,1),'-','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.6,0.1,0.9];  plot(x,INFO(:,ii,8,1),'-','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.6,0.3,0.0];  plot(x,INFO(:,ii,9,1),'-','Color',cl, 'MarkerFaceColor',cl);
    
    ylabel('Information (bit)');
    set(gca,'Ytick',0:0.2:1,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'}, 'TickLength',[0.01,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.01,0]);
    axis([0 27.5 0 1]);  box off;  % for plotting measured noise
    set(gca,'FontSize',10);
end

set(gcf,'Color',[1,1,1]);

% % save figure
% export_fig('PerfInfo.pdf')

end