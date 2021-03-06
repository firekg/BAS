function plot_info_progress_v2(INFOA, parS)

[parM, INFOAC] = combine_parS_INFOA(INFOA, parS);

figure(5);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

pd=0.04; %distance between plots
pws=0.23; %plot width
phs=0.23; %plot height
lm=0.08; %left margin
ls1=lm; %1st column left start position
ts1=0.75; %1st row top start position

nss = 5+6;  %number of session
nboot = 50;  %number of bootstraps

%INFOAC has 5 layers: rev num; parti; session; moments; boot;
INFO = nanmean(INFOAC,5);
x = 1:25;
curveM = zeros(25,nss,nboot);
for ii=1:3 %column
for jj=1  %row
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts1-(jj-1)*phs-(jj-1)*pd;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
        
    %parM has 3 layers: iboot; npar; parti
    for iboot = 1:nboot
        par = parM(iboot,:,ii);
        curveM(:,:,iboot) = showfit(par, nss);
    end
    curves = mean(curveM,3);
    
    for ss = 1:nss
        % grayscale color map
        %cl = (0.08*ss)*ones(1,3); %(0.08*ss + 0.12)*ones(1,3);
        % rainbow color map
        cl(1) =  (nss-ss)/(nss-1); %red
        if (ss <= nss/2) %green
            cl(2) = (ss-1)/(nss-1)*2;
        else
            cl(2) = (nss-ss)/(nss-1)*2;
        end
        cl(3) = (ss-1)/(nss-1); %blue
        % plot
        plot(x,INFO(:,ii,ss,1),'-','Color',cl);  hold on;
        plot(x,curves(:,ss),':','Color',cl);
    end
    
    ylabel('Information (bit)');
    set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1'}, 'TickLength',[0.02,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.02,0]);
    axis([0 27.5 0 0.8]);  box off;  
    set(gca,'FontSize',10);
    
%     legContent = cast(reshape(repmat(1:11,2,1),22,1), 'char');
%     legend(legContent);
end
end

for ii=1:3 %column
for jj=2  %row
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts1-(jj-1)*phs-(jj-1)*pd;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);

    for ss = 1:nss
        cl(1) =  (nss-ss)/(nss-1); %red
        if (ss <= nss/2) %green
            cl(2) = (ss-1)/(nss-1)*2;
        else
            cl(2) = (nss-ss)/(nss-1)*2;
        end
        cl(3) = (ss-1)/(nss-1); %blue
        % plot
        xx = ss+randn(nboot, 1)*0.05;
        base = parM(:, nss+1, ii); %suppose to be largest
        yy = base./parM(:, ss+1, ii);
        %plot(xx, yy, '.', 'Color', cl, 'MarkerSize', 1e-3); hold on;
        plot(ss, mean(yy), 'o', 'Color', cl); hold on;
        errorbar(ss, mean(yy), mad(yy)*1.4826, 'Color', cl);
    end
    plot([5.5, 5.5], [0, 2], 'k:');
    plot([0, 100], [1, 1], 'k:');
    
    axis([0 nss+1 0.5 1.5]);
    box off;  
    ylabel('Relative efficiency');
    %set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','','','','','0.5','','','','','1'}, 'TickLength',[0.02,0]);    
    xlabel('Session number');
    set(gca,'Xtick',0:1:11,'XTickLabel',{'','1','2','3','4','5','1','2','3','4','5','6'}, 'TickLength',[0.02,0]);
    set(gca,'FontSize',10);
    
    yy = parM(:, 2:nss+1, ii);
    day1 = mean(mean(yy(:,1:5),1),2);
    day2n3 = mean(mean(yy(:,6:11),1),2);
    fss3n4 = mean(mean(yy(:,8:9),1),2);
    fss5n6 = mean(mean(yy(:,10:11),1),2);
    fprintf('parti %d: day 1 scale = %f; day 2&3 scale = %f; free-scan 3&4 = %f; free-scan 5&6 = %f;\n',...
        ii, day1, day2n3, fss3n4, fss5n6);
    fprintf('day 2&3 / 1 = %f; free-scan 5&6 / 3&4 = %f;\n',...
        day2n3/day1, fss5n6/fss3n4);

end
end

set(gcf,'Color',[1,1,1]);
% export_fig('ALprogress.pdf')

end

%% quadratic fitfunc
function curves = showfit(par, nss)

x = (1:25)';
curves = zeros(25,nss);
for ss = 1:nss
    curves(:,ss) = wblcdf(x, par(1+ss), par(1));
end

end

%% combine things 
function [parMC, INFOAC] = combine_parS_INFOA(INFOA, parS)
%INFOActive has 5 layers: rev num; parti; session; moments; iboot;
%parM has 3 layers: iboot; npar; parti

nboot = 5;
nrun = 10;
nss = 5+6;  %number of session
npar = 1+nss;  %number of fit parameters
nparti = 3;

parMC = zeros(nboot*nrun, npar, nparti);
INFOAC = zeros(25, nparti, nss, 2, nboot*nrun);

% start combining data
for irun = 1:nrun
    parMC((irun-1)*nboot+1:irun*nboot ,:,:) = parS(irun).par;
    INFOAC(:,:,:,:,(irun-1)*nboot+1:irun*nboot) = INFOA(irun).info;
end

end


