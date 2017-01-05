function [varargout] = info_progress_v2(DALwPix)

ns.ss = 11; %5+6  %number of session
ns.par = 1+ns.ss;  %number of fit parameters
ns.boot = 5; %5  %number of bootstraps
ns.samp = 200; %20  %number of perception-integration samples
ns.parti = 3;

INFOC = prop_progress(DALwPix, ns);
%INFOC = info_prop(PROP, ns); % can delete
INFOActive = nanmean(INFOC,6);  %INFOActive has 5 layers: rev num; parti; session; moments; iboot;

%parM has 3 layers: iboot; npar; parti
parM = zeros(ns.boot, ns.par, ns.parti);
for parti= 1:ns.parti
    for iboot = 1:ns.boot
        parInit = 15*ones(1, ns.par); % initialize wblcdf scale parameter at 15
        parInit(1) = 1; % initialize wblcdf power parameter at 1
        INFOboot = INFOActive(:,parti,:,:,iboot);
        myfun = @(par)er2_info(par, INFOboot, ns.ss);
        options = optimset('MaxFunEvals', 50000, 'MaxIter', 50000);
        [parM(iboot,:,parti),fval,exitflag] = fminsearch(myfun, parInit, options);
        fprintf('fval=%f; flag=%f.\n', fval, exitflag);
    end
end
% varargout = {INFOActive, parM};
varargout = {INFOActive, parM};

end
%%
function [INFO] = prop_progress(DALwPix, ns)

    INFO = nan(25, ns.parti, ns.ss, 2, ns.boot, ns.samp);

    %pv= [  lsigo,    phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, window, prxpa, prCommon, shift];
    pv =[log(1.0),   nan,     nan,   nan,       log(0),        nan,     0,  1/2,      1,   1/3,      1/3,    12];
    %pw=[lbm,  km,  nm, bias];
    pw=[   0, nan, nan,    0]; % none of the last 3 are used here
    
    % make probe points z
    npix=770;
    %ImageSize=20;
    n=110;
    xx = linspace(1,npix,n);
    yy = linspace(1,npix,n);
    [x1,x2]=meshgrid(xx,yy);
    z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
    %dxx = xx(2)-xx(1);
        
    for parti = 1:ns.parti %3
        if(parti==1)
            pv(1) = log(0.5);
            pv(6) = 0.34901995;
            pv(7) = 0.04471412;
            pv(12) = 16;
        elseif(parti==2)  
            pv(1) = log(0.5);
            pv(6) = 0.637376715;
            pv(7) = 0.120818875;
            pv(12) = 17;
        elseif(parti==3)  
            pv(1) = log(0.3);
            pv(6) = 0.374908933;
            pv(7) = 0.101003052;
            pv(12) = 15;
        end
        for ss = 1:ns.ss
        D = DALwPix(parti, ss);
            for iboot = 1:ns.boot
                %bootstrap with replacement each trial
                nTr = D.Trials;
                ind = ceil(rand(nTr,1)*nTr);
                AFboot = D;
                AFboot.RevealPosX = D.RevealPosX(ind,:);
                AFboot.RevealPosY = D.RevealPosY(ind,:);
                AFboot.RevealZ = D.RevealZ(ind,:);
                AFboot.ImageID = D.ImageID(ind,:);
                AFboot.MaxRevealingTrial = D.MaxRevealingTrial(ind,:);
                Gmap = zeros(nTr,25,3);
                for samp = 1:ns.samp       
                    PROP = a_Get_BALDscoreProp(AFboot,[],pv,pw,randn(1000,25),Gmap,z,2);
                    % needed below to save memory explosion
                    m = PROP.mla;  
                    M = -log(0.5) +m.*log(m) +(1-m).*log(1-m);  
                    M = -M/log(0.5);
                    INFO(:, parti, ss, 1, iboot, samp) = nanmean(M);  
                    INFO(:, parti, ss, 2, iboot, samp) = sqrt(nanvar(M)./sum(~isnan(M)));
                end
            end    
        fprintf('Done parti %d, session %d.\n', parti, ss);
        end
    end

end

%% Compute information gain -- can delete
% function [INFO] = info_prop(PROP, ns)
% % PROP has 4 layers: iboot; samp; session; parti
% % INFO has 6 layers: revealing number (1:25); parti; session;
% %                    moments x2 (mean,SD); nboot; samp
% INFO = nan(25, ns.parti, ns.ss, 2, ns.boot, ns.samp);
% for iboot = 1:ns.boot
%     for samp = 1:ns.samp
%         for ss = 1:ns.ss
%             for parti = 1:ns.parti
%                 m = PROP(parti, ss, iboot, samp).mla;  
%                 M = -log(0.5) +m.*log(m) +(1-m).*log(1-m);  
%                 M = -M/log(0.5);
%                 INFO(:, parti, ss, 1, iboot, samp) = nanmean(M);  
%                 INFO(:, parti, ss, 2, iboot, samp) = sqrt(nanvar(M)./sum(~isnan(M)));
%             end
%         end
%     end
% end
% end

%% 
function er2 = er2_info(par, INFOiboot, nss)

x = (1:25)';
er2 = 0;
%INFOiboot has 4 layers: rev num; dummy; session; moments;
for ss = 1:nss
   er2 = er2 + norm( wblcdf(x, par(1+ss), par(1)) - INFOiboot(:,1,ss,1) ); 
end

end






