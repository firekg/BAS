figure(1);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[0,0,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [0,0,20,20]);

% subplot('Position',[0,0.7,0.3,0.3]);
% plot(1:100,1:100,'ro');
% set(gca,'Units','centimeters');
% set(gca,'Position',[1,1,8,4]);
% % legend('a','Location','NorthEastOutside');
% set(gca,'Xtick',[1,5,25,97],'XTickLabel',{'a','b','c','d'});
% xlabel('Letters'); ylabel('Numbers');
% ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);
% xlabh = get(gca,'XLabel'); set(xlabh,'Position',get(xlabh,'Position') + [0, 3, 0]);

% ##########################################################################
d=0.03;
load('SCdata.mat');
%
subplot('Position',[0+2*d,0.7,0.15-d/2,0.3-d]);
    errorbar(x,PerfPLB,sdPLB,'k*-');  hold on;
    errorbar(x,PerfPLaB,sdPLaB,'+-','Color',[0.5,0.5,0.5]);
    errorbar(x,PerfPLR,sdPLR,'bs-');
    errorbar(x,PerfAL,sdAL,'ro-');
set(gca,'FontSize',9);
% xlabel('N revealing');
ylabel('Performance');
% ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);
axis([0 30 0.4 1]);
set(gca,'Xtick',0:5:30,'XTickLabel',{'','','10','','20','',''});
title('Subject')

subplot('Position',[0.15+1.5*d,0.7,0.15-d/2,0.3-d]);
    errorbar(x,mPerfPLB,msdPLB,'k*-');  hold on;
    errorbar(x,mPerfPLaB,msdPLaB,'+-','Color',[0.5,0.5,0.5]);
    errorbar(x,mPerfPLR,msdPLR,'bs-');
    errorbar(x,mPerfAL,msdAL,'ro-');
set(gca,'FontSize',9);
% xlabel('N revealing');
ylabel('');
set(gca,'YTickLabel','');
axis([0 30 0.4 1]);
set(gca,'Xtick',0:5:30,'XTickLabel',{'','','10','','20','',''});
title('Model')

subplot('Position',[0.3+3*d,0.7,0.3-3*d,0.3-d]);
    Bx=1:25;
%     errorbar(Bx, SCPLBptile_m, sqrt(SCPLBptile_v./SCPLBptile_n), 'k-*'); hold on;
%     errorbar(Bx, SCPLaBptile_m, sqrt(SCPLaBptile_v./SCPLaBptile_n), '-+', 'Color',[0.5,0.5,0.5]);
    errorbar(Bx, SCPLRptile_m, sqrt(SCPLRptile_v./SCPLRptile_n), 'b-o', 'MarkerFaceColor', [0,0,1]); hold on;
    errorbar(Bx, SCPLBsoftptile_m, sqrt(SCPLBsoftptile_v./SCPLBsoftptile_n), 'g-o', 'MarkerFaceColor', [0,1,0]);
    errorbar(Bx, SCALptile_m, sqrt(SCALptile_v./SCALptile_n), 'r-o', 'MarkerFaceColor', [1,0,0]);  hold off;
set(gca,'FontSize',9);
% xlabel('N Revealing');
ylabel('BAS score percentile');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [1.5, 0, 0]);
axis([0 30 40 100]);  %axis square;  %grid on;  
set(gca,'Xtick',0:5:30,'XTickLabel',{'0','','10','','20','','30'});
set(gca,'Ytick',40:20:100);
% legend('BALD', 'anti-BALD', 'random', 'AL', 'AL model', 'Location','SouthWest'); legend('boxoff');

%
subplot('Position',[0.6+d,0.85-d/3,0.12,0.12]);
image((dhrmpa-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('patchy')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);
ylabel('Subject');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [30, 0, 0]);

subplot('Position',[0.72+d,0.85-d/3,0.12,0.12]);
image((dhrmsh-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('hor str')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);

subplot('Position',[0.84+d,0.85-d/3,0.12,0.12]);
image((dhrmsv-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('ver str')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);

%
subplot('Position',[0.6+d,0.7+d/3,0.12,0.12]);
image((dmrmpa-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
ylabel('Model');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [30, 0, 0]);

subplot('Position',[0.72+d,0.7+d/3,0.12,0.12]);
image((dmrmsh-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

subplot('Position',[0.84+d,0.7+d/3,0.12,0.12]);
image((dmrmsv-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
%
clear

% ##########################################################################
d=0.03;
load('EPdata.mat');
%
subplot('Position',[0+2*d,0.4,0.15-d/2,0.3-d]);
    errorbar(x,PerfPLB,sdPLB,'k*-');  hold on;
    errorbar(x,PerfPLaB,sdPLaB,'+-','Color',[0.5,0.5,0.5]);
    errorbar(x,PerfPLR,sdPLR,'bs-');
    errorbar(x,PerfAL,sdAL,'ro-');
set(gca,'FontSize',9);
% xlabel('N revealing');
ylabel('Performance');
% ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);
axis([0 30 0.4 1]);
set(gca,'Xtick',0:5:30,'XTickLabel',{'','','10','','20','',''});
% title('Subject')

subplot('Position',[0.15+1.5*d,0.4,0.15-d/2,0.3-d]);
    errorbar(x,mPerfPLB,msdPLB,'k*-');  hold on;
    errorbar(x,mPerfPLaB,msdPLaB,'+-','Color',[0.5,0.5,0.5]);
    errorbar(x,mPerfPLR,msdPLR,'bs-');
    errorbar(x,mPerfAL,msdAL,'ro-');
set(gca,'FontSize',9);
% xlabel('N revealing');
ylabel('');
set(gca,'YTickLabel','');
axis([0 30 0.4 1]);
set(gca,'Xtick',0:5:30,'XTickLabel',{'','','10','','20','',''});
% title('Model')

subplot('Position',[0.3+3*d,0.4,0.3-3*d,0.3-d]);
    Bx=1:25;
%     errorbar(Bx, PLBptile_m, sqrt(PLBptile_v./PLBptile_n), 'k-*'); hold on;
%     errorbar(Bx, PLaBptile_m, sqrt(PLaBptile_v./PLaBptile_n), '-+', 'Color',[0.5,0.5,0.5]);
    errorbar(Bx, PLRptile_m, sqrt(PLRptile_v./PLRptile_n), 'b-o', 'MarkerFaceColor', [0,0,1]);  hold on;
    errorbar(Bx, PLBsoftptile_m, sqrt(PLBsoftptile_v./PLBsoftptile_n), 'g-o', 'MarkerFaceColor', [0,1,0]);
    errorbar(Bx, ALptile_m, sqrt(ALptile_v./ALptile_n), 'r-o', 'MarkerFaceColor', [1,0,0]);  hold off;
set(gca,'FontSize',9);
% xlabel('N Revealing');
ylabel('BAS score percentile');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [1.5, 0, 0]);
axis([0 30 40 100]);  %axis square;  %grid on;  
set(gca,'Xtick',0:5:30,'XTickLabel',{'0','','10','','20','','30'});
set(gca,'Ytick',40:20:100);
% legend('BALD', 'anti-BALD', 'random', 'AL', 'AL model', 'Location','SouthWest'); legend('boxoff');

%
subplot('Position',[0.6+d,0.55-d/3,0.12,0.12]);
image((dhrmpa-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('patchy')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);
ylabel('Subject');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [30, 0, 0]);

subplot('Position',[0.72+d,0.55-d/3,0.12,0.12]);
image((dhrmsh-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('hor str')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);

subplot('Position',[0.84+d,0.55-d/3,0.12,0.12]);
image((dhrmsv-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('ver str')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);

%
subplot('Position',[0.6+d,0.4+d/3,0.12,0.12]);
image((dmrmpa-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
ylabel('Model');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [30, 0, 0]);

subplot('Position',[0.72+d,0.4+d/3,0.12,0.12]);
image((dmrmsh-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

subplot('Position',[0.84+d,0.4+d/3,0.12,0.12]);
image((dmrmsv-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
%
clear

% ########################################################################## 
d=0.03;
load('BZdata.mat');
%
subplot('Position',[0+2*d,0.1,0.15-d/2,0.3-d]);
    errorbar(x,PerfPLB,sdPLB,'k*-');  hold on;
    errorbar(x,PerfPLaB,sdPLaB,'+-','Color',[0.5,0.5,0.5]);
    errorbar(x,PerfPLR,sdPLR,'bs-');
    errorbar(x,PerfAL,sdAL,'ro-');
set(gca,'FontSize',9);
xlabel('N revealing');
ylabel('Performance');
% ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);
axis([0 30 0.4 1]);
set(gca,'Xtick',0:5:30,'XTickLabel',{'','','10','','20','',''});
% title('Subject')

subplot('Position',[0.15+1.5*d,0.1,0.15-d/2,0.3-d]);
    errorbar(x,mPerfPLB,msdPLB,'k*-');  hold on;
    errorbar(x,mPerfPLaB,msdPLaB,'+-','Color',[0.5,0.5,0.5]);
    errorbar(x,mPerfPLR,msdPLR,'bs-');
    errorbar(x,mPerfAL,msdAL,'ro-');
set(gca,'FontSize',9);
xlabel('N revealing');
ylabel('');
set(gca,'YTickLabel','');
axis([0 30 0.4 1]);
set(gca,'Xtick',0:5:30,'XTickLabel',{'','','10','','20','',''});
% title('Model')

subplot('Position',[0.3+3*d,0.1,0.3-3*d,0.3-d]);
    Bx=1:25;
%     errorbar(Bx, PLBptile_m, sqrt(PLBptile_v./PLBptile_n), 'k-*'); hold on;
%     errorbar(Bx, PLaBptile_m, sqrt(PLaBptile_v./PLaBptile_n), '-+', 'Color',[0.5,0.5,0.5]);
    errorbar(Bx, PLRptile_m, sqrt(PLRptile_v./PLRptile_n), 'b-o', 'MarkerFaceColor', [0,0,1]);  hold on;
    errorbar(Bx, PLBsoftptile_m, sqrt(PLBsoftptile_v./PLBsoftptile_n), 'g-o', 'MarkerFaceColor', [0,1,0]);
    errorbar(Bx, ALptile_m, sqrt(ALptile_v./ALptile_n), 'r-o', 'MarkerFaceColor', [1,0,0]);  hold off;
set(gca,'FontSize',9);
xlabel('Revealing number');
ylabel('BAS score percentile');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [1.5, 0, 0]);
axis([0 30 40 100]);  %axis square;  %grid on;  
set(gca,'Xtick',0:5:30,'XTickLabel',{'0','','10','','20','','30'});
set(gca,'Ytick',40:20:100);
% legend('BALD', 'anti-BALD', 'random', 'AL', 'AL model', 'Location','SouthWest'); legend('boxoff');

%
subplot('Position',[0.6+d,0.25-d/3,0.12,0.12]);
image((dhrmpa-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('patchy')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);
ylabel('Subject');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [30, 0, 0]);

subplot('Position',[0.72+d,0.25-d/3,0.12,0.12]);
image((dhrmsh-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('hor str')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);

subplot('Position',[0.84+d,0.25-d/3,0.12,0.12]);
image((dhrmsv-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title('ver str')
tlabh = get(gca,'Title'); set(tlabh,'Position',get(tlabh,'Position') + [0, 40, 0]);

%
subplot('Position',[0.6+d,0.1+d/3,0.12,0.12]);
image((dmrmpa-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
ylabel('Model');
ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [30, 0, 0]);

subplot('Position',[0.72+d,0.1+d/3,0.12,0.12]);
image((dmrmsh-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

subplot('Position',[0.84+d,0.1+d/3,0.12,0.12]);
image((dmrmsv-vmin)/dv*(nc-1)+1);
colormap(hot);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
%
% clear

% saveas(gcf,'3subjects','pdf')

