function plot_info_progress()

[parM, INFOActive] = combine_parM_INFOActive();
%[parMA] = parM_INFOasym(INFOActive);

figure(5);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
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

%INFOActive has 5 layers: rev num; parti; session; moments; iboot;
INFO = nanmean(INFOActive,5);
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
        cl = (0.08*ss)*ones(1,3); %(0.08*ss + 0.12)*ones(1,3);
        plot(x,INFO(:,ii,ss,1),'-','Color',cl);  hold on;
        plot(x,curves(:,ss),':','Color',cl);
    end
    
    ylabel('Information (bit)');
    set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','','','','','1'}, 'TickLength',[0.02,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.02,0]);
    axis([0 27.5 0 0.31]);  box off;  
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

    %parM has 3 layers: iboot; npar; parti
    xFF = repmat(1:5, nboot, 1);
    yFF = parM(:, 3:7, ii)./1e-3;    
    [rho, pval] = corr(xFF(:), yFF(:), 'type', 'Pearson');
    fprintf('Parti %d Familiarize:  r=%f; p=%f.\n', ii, rho, pval);
    xFS = repmat(6:nss, nboot, 1);
    yFS = parM(:, 8:end, ii)./1e-3;    
    [rho, pval] = corr(xFS(:), yFS(:), 'type', 'Pearson');
    fprintf('Parti %d Free-scan:    r=%f; p=%f.\n', ii, rho, pval);
%     xFS = repmat(9:nss, nboot, 1);
%     yFS = parM(:, 11:end, ii)./1e-3;    
%     [rho, pval] = corr(xFS(:), yFS(:), 'type', 'Pearson');
%     fprintf('Parti %d last 3 Free-scan:    r=%f; p=%f.\n', ii, rho, pval);
    xFA = repmat(1:nss, nboot, 1);
    yFA = parM(:, 3:end, ii)./1e-3;    
    [rho, pval] = corr(xFA(:), yFA(:), 'type', 'Pearson');
    fprintf('Parti %d Overall:      r=%f; p=%f.\n', ii, rho, pval);

    % ratio between 1st session and the average of last 6 sessions
    yFF1 = parM(:, 3, ii)./1e-3;
    ss1 = mean(yFF1);
    yFS = parM(:, 8:end, ii)./1e-3;
    ssFS = mean(mean(yFS,2));
    FSoverFF1 = mean(yFS,2)./yFF1;
    FSoverFF1 = sort(FSoverFF1);
    fprintf('Parti %d ssFS/ss1 =    %f; \n', ii, ssFS/ss1 );
%     [h,p] = ttest2(yFF1(:), yFS(:));
%     fprintf('Parti %d: ttest h=%f; p=%f; \n', ii, h, p);
%     fprintf('Parti %d: ratio 95 CI: %f--%f; \n\n', ii, FSoverFF1(1), FSoverFF1(end));

    % ratio between 5&6 sessions and 3&4 sessions
    yFS34 = parM(:, 10:11, ii)./1e-3;
    ss34 = mean(mean(yFS34,2));
    yFS56 = parM(:, 12:13, ii)./1e-3;
    ss56 = mean(mean(yFS,2));
    fprintf('Parti %d ss56/ss34 =    %f; \n', ii, ss56/ss34 );


%     %parMA has 3 layers: npar; parti; iboot (different from parM!)
%     % ratio between sessions 10-11/8-9
%     ssPar = parMA(:, ii, :)./1e-3;
%     ssPar = mean(ssPar,3);
%     fprintf('Parti %d ss10-11/8-9 = %f; \n', ii, ssPar(4)/ssPar(3));

    fprintf('\n');

    xx = repmat(1:nss, nboot, 1) + randn(nboot, nss)*0.05;
    yy = parM(:, 3:end, ii)./1e-3;
    plot(xx(:), yy(:), 'k.', 'MarkerSize', 1e-3); hold on;
    plot([5.5, 5.5], [0, 1], 'k:');
    
    axis([0 nss+1 0 1]);  box off;  
    ylabel('Scale (a.u.)');
    set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','','','','','0.5','','','','','1'}, 'TickLength',[0.02,0]);    
    xlabel('Session number');
    set(gca,'Xtick',0:1:11,'XTickLabel',{'','1','2','3','4','5','1','2','3','4','5','6'}, 'TickLength',[0.02,0]);
    set(gca,'FontSize',10);
end
end

set(gcf,'Color',[1,1,1]);


end

%% quadratic fitfunc
function curves = showfit(par, nss)

x = (1:25)';
curves = zeros(25,nss);
for ss = 1:nss
    curves(:,ss) = par(2+ss)*(x.^2 + par(1)*x + par(2));
end

end

%% combine things 
function [parMC, INFOActiveC]=combine_parM_INFOActive()
%INFOActive has 5 layers: rev num; parti; session; moments; iboot;
%parM has 3 layers: iboot; npar; parti

nboot = 5;
nrun = 10;
nss = 5+6;  %number of session
npar = 2+nss;  %number of fit parameters
nparti = 3;

parMC = zeros(nboot*nrun, npar, nparti);
INFOActiveC = zeros(25, nparti, nss, 2, nboot*nrun);

% start combining data
for irun = 1:nrun
    load(['InfoProgress_data/parM',num2str(irun),'.mat']);
    parMC((irun-1)*nboot+1:irun*nboot ,:,:) = parM;
    INFOActiveC(:,:,:,:,(irun-1)*nboot+1:irun*nboot) = INFOActive;
end

end

%% combine things 2
function [parM]=parM_INFOasym(INFOActiveC)
%INFOActive has 5 layers: rev num; parti; session; moments; iboot;
%parM has 3 layers: npar; parti; iboot
nboot = 50;
nss = 2;
npar = 2+nss;  %number of fit parameters
nparti = 3;

INFOasym(:,:,1,:,:) = mean(INFOActiveC(:,:,8:9,:,:),3);  % 8-9 combined
INFOasym(:,:,2,:,:) = mean(INFOActiveC(:,:,10:11,:,:),3);  % 10-11 combined

%parM has 3 layers: iboot; npar; parti (copied from info_progress)
parM = zeros(npar, nparti, nboot);
for parti= 1:nparti
    for iboot = 1:nboot
        par0 = zeros(npar,1);
        INFOiboot = INFOasym(:,parti,:,:,iboot);
        myfun = @(par)er2_info(par, INFOiboot, nss);
        options = optimset('MaxIter', 1000);
        [parM(:,parti,iboot),fval,exitflag] = fminsearch(myfun, par0, options);
        %fprintf('fval=%f; flag=%f.\n', fval, exitflag);
    end
end

end

%% copied from info_progress
function er2 = er2_info(par, INFOiboot, nss)

x = (1:25)';
er2 = 0;
%INFOiboot has 4 layers: rev num; dummy; session; moments;
for ss = 1:nss
   er2 = er2 + norm( par(2+ss)*(x.^2 + par(1)*x + par(2)) - INFOiboot(:,1,ss,1)); 
end

end

