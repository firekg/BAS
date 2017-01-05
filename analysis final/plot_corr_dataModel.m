function plot_corr_dataModel(CORR)
%%
figure(4);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

pd=0.04; %distance between plots
pws=0.18; %plot width
phs=0.15; %plot height
lm=0.08; %left margin
ls1=lm; %1st column left start position
ts1=0.75; %1st row top start position
x=1:25;
ticklen = 0.02;

% CORR has 5 layers:
% moments x3: (mean, upper 95%, lower 95%)
% revealing number 25x: 1:25
% parti-cond x2: within parti, across parti
% patt-cond x2: same patt, diff parti
% participants x3: SC, EP, BZ

for parti=1:3 
for icond=1:2 %parti-cond
    ls = ls1+(parti-1)*pws+(parti-1)*pd;
    ts = ts1-(icond-1)*phs-(icond-1)*pd;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
    plot([0,30],[0,0],'-','Color',[0.0,0.0,0.0]); hold on;   %line on 0
    % same patt
    cl=[1.0,0.5,0.0];
    plot(x, CORR(1,:,icond,1,parti),'-','Color',cl, 'MarkerFaceColor',cl);
    a_ErrorShade(x, CORR(2,:,icond,1,parti), CORR(3,:,icond,1,parti),cl,1);
    % diff patt
    cl=[0.5,0.0,1.0];  
    plot(x, CORR(1,:,icond,2,parti),'-','Color',cl, 'MarkerFaceColor',cl);
    a_ErrorShade(x, CORR(2,:,icond,2,parti), CORR(3,:,icond,2,parti),cl,1);

    ylabel('Correlation');
    set(gca,'Ytick',-1:0.5:1,'YTickLabel',{'-1','','0','','1'}, 'TickLength',[ticklen,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[ticklen,0]);
    axis([0 25 -1.0 1.0]);  box off;  
    set(gca,'FontSize',10);
    
    if (icond==1)
        fprintf('parti %d with self, same patt: %f(%f-%f).\n',...
            parti, CORR(1,25,icond,1,parti), CORR(3,25,icond,1,parti), CORR(2,25,icond,1,parti));
        fprintf('parti %d with self, diff patt: %f(%f-%f).\n',...
           parti, CORR(1,25,icond,2,parti), CORR(3,25,icond,2,parti), CORR(2,25,icond,2,parti));
    elseif(icond==2)
        fprintf('parti %d with others, same patt: %f(%f-%f).\n',...
            parti, CORR(1,25,icond,1,parti), CORR(3,25,icond,1,parti), CORR(2,25,icond,1,parti));
        fprintf('parti %d with others, diff patt: %f(%f-%f).\n',...
           parti, CORR(1,25,icond,2,parti), CORR(3,25,icond,2,parti), CORR(2,25,icond,2,parti));
    end
end
end

set(gcf,'Color',[1,1,1]);

% save figure
% export_fig('CorrModel.pdf')


end