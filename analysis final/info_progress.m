function [varargout]=info_progress(DActive, INFO)

if (nargin==1) %the second argin is not used
    nss = 5+6;  %number of session
    npar = 2+nss;  %number of fit parameters
    nboot = 5;  %number of bootstraps
    nsamp = 20;  %number of perception-integration samples
    nparti = 3;
    
    [PROP] = info_AF(DActive, nss, nboot, nsamp);
    [INFOC] = info_prop(PROP, nss, nboot, nsamp);
    INFOActive = nanmean(INFOC,6);  %INFOActive has 5 layers: rev num; parti; session; moments; iboot;

    %parM has 3 layers: iboot; npar; parti
    parM = zeros(nboot, npar, nparti);
    for parti= 1:nparti
        for iboot = 1:nboot
            par0 = zeros(1, npar);
            INFOiboot = INFOActive(:,parti,:,:,iboot);
            myfun = @(par)er2_info(par, INFOiboot, nss);
            options = optimset('MaxIter', 1000);
            [parM(iboot,:,parti),fval,exitflag] = fminsearch(myfun, par0, options);
            fprintf('fval=%f; flag=%f.\n', fval, exitflag);
            if (0) % plot to check
                curves = showfit(parM(iboot,:));
                for ss = 1:nss
                    plot((1:25)', INFOiboot(:,1,ss,1), '-'); hold on;
                    plot((1:25)', curves(:,ss), ':'); hold off;        
                    pause
                end
            end
        end
    end
    varargout = {INFOActive, parM};
elseif (nargin==2)  %the first argin is not used
    ncond = 4;  %number of conditions: AL, PLR, PLB, nBAS, p-indep & order-dep, p-dep & order-indep, p-dep & order-dep
    npar = 2+ncond;  %number of fit parameters
    nboot = 200;
    %parM has 3 layers: npar; parti, iboot
    parM = zeros(npar, 4, nboot);
    for iboot = 1:nboot
        for parti= 1:4
            par0 = zeros(1, npar);
            % INFO has 5 layers: rev num; parti x4; cond x6 (AL,PLR,PLB,...); moments x2 (mean,SD);
            INFOi = INFO(:, parti, :, 1, iboot);
            myfun = @(par)er2_info(par, INFOi, ncond);
            options = optimset('MaxIter', 1000);
            [parM(:, parti, iboot), fval, exitflag] = fminsearch(myfun, par0, options);
            if (0) % plot to check
                curves = showfit(parM(iboot,:));
                for ss = 1:nss
                    plot((1:25)', INFOi(:,1,cond,1), '-'); hold on;
                    plot((1:25)', curves(:,cond), ':'); hold off;        
                    pause
                end
            end
        end
    end
    [ratioM] = info_ratio(parM);
    varargout = {parM};    
end

end
%% compute info
function [PROP]=info_AF(DActive, nss, nboot, nsamp)
% DActive has 2 layers: parti x3; seesions
%PROP has 4 layers: iboot; samp; session; parti

for parti = 1:3
    if (parti==1)       sigo = 1.1;
    elseif (parti==2)   sigo = 1.2;
    elseif (parti==3)   sigo = 0.9;
    end
    for ss = 1:nss;
        fprintf('\nSession %d/%d.\n', ss, nss);
        D = DActive(parti, ss);
        fprintf(' Bootstrap sample (%d): ', nboot);
        PROP(:,:,ss,parti) = prop_AF(D, sigo, nboot, nsamp);  %change sigo
    end
    fprintf('\nParti %d complete.\n', parti);
end

end
%% compute PROP of active familiarization trials
function [PROP] = prop_AF(DA, sigo, nboot, nsamp)
% PROP has 2 layer: iboot, samp;
    Gmap=1;
    %pv= [ lsigo,  phit,   lphis, phis0,   lbx,    lbd, lapse, prpa, w];
    pv =[log(1),   nan,     nan,   nan,    nan,    nan,     0,  0.5, 1];
    %pw=[lbm,  km,  nm, bias];
    pw=[   0, nan, nan,    0]; % none of the last 3 are used here
    % make probe points z
    npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
    [x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
    dxx = xx(2)-xx(1);
    
    DIM = a_DIM(DA);    
    
    for iboot = 1:nboot
        %bootstrap with replacement each trial
        nTr = DA.Trials;
        ind = ceil(rand(nTr,1)*nTr);
        AFboot = DA;
        AFboot.RevealPosX = DA.RevealPosX(ind,:);
        AFboot.RevealPosY = DA.RevealPosY(ind,:);
        AFboot.ImageID = DA.ImageID(ind,:);
        AFboot.MaxRevealingTrial = DA.MaxRevealingTrial(ind,:);
    
        DIMboot = DIM(ind);
        pv(1) = log(sigo);
    
        for samp = 1:nsamp        
            PROP(iboot, samp) = a_Get_BALDscoreProp(AFboot,DIMboot,pv,pw,randn(1000,25),Gmap,z,2);
        end
        fprintf('%d, ', iboot);
    end
end

%% Compute information gain
function [INFO] = info_prop(PROP, nss, nboot, nsamp)
% PROP has 4 layers: iboot; samp; session; parti
% INFO has 6 layers: revealing number (1:25); parti; session;
%                    moments x2 (mean,SD); nboot; samp
INFO = nan(25, 3, nss, 2, nboot, nsamp);
for iboot = 1:nboot
    for samp = 1:nsamp
        for ss = 1:nss
            for parti = 1:3
                m=PROP(iboot, samp, ss, parti).mla;  
                M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  
                M= -M/log(0.5);
                INFO(:, parti, ss, 1, iboot, samp) = nanmean(M);  
                INFO(:, parti, ss, 2, iboot, samp) = sqrt(nanvar(M)./sum(~isnan(M)));
            end
        end
    end
end
end
%% 
function er2 = er2_info(par, INFOiboot, nss)

x = (1:25)';
er2 = 0;
%INFOiboot has 4 layers: rev num; dummy; session; moments;
for ss = 1:nss
   er2 = er2 + norm( par(2+ss)*(x.^2 + par(1)*x + par(2)) - INFOiboot(:,1,ss,1)); 
end

end
%% quadratic fitfunc
function curves = showfit(par, nss)

x = (1:25)';
curves = zeros(25,nss);
for ss = 1:nss
    curves(:,ss) = par(2+ss)*(x.^2 + par(1)*x + par(2));
end

end

%% given parMC, get ratios
function [ratioM] = info_ratio(parMC)

% ratioM has two layers: mean & UB & LB; parti; cond (AL/PLR, PLB/AL, BAS/AL);
nRatio = 6;
ratioM = zeros(3, 4, nRatio);
for iparti = 1:4
    
    iRatio = 1; %cond AL/PLR
    temp = parMC(3,iparti,:)./parMC(4,iparti,:);
    temp = sort(temp);
    ratioM(1,iparti,iRatio) = mean(temp);
    ratioM(2,iparti,iRatio) = temp(190);
    ratioM(3,iparti,iRatio) = temp(10);
    fprintf('parti %d: AL/rand=%f; Upper=%f; Lower=%f\n',...
        iparti, ratioM(1,iparti,iRatio), ratioM(2,iparti,iRatio), ratioM(3,iparti,iRatio));
    
    iRatio = 2; %cond PLB/AL
    temp = parMC(5,iparti,:)./parMC(3,iparti,:);
    temp = sort(temp);
    ratioM(1,iparti,iRatio) = mean(temp);
    ratioM(2,iparti,iRatio) = temp(190);
    ratioM(3,iparti,iRatio) = temp(10);
    fprintf('parti %d: optimal/AL=%f; Upper=%f; Lower=%f\n',...
        iparti, ratioM(1,iparti,iRatio), ratioM(2,iparti,iRatio), ratioM(3,iparti,iRatio));
    
    iRatio = 3; %cond BAS/AL
    temp = parMC(6,iparti,:)./parMC(3,iparti,:);
    temp = sort(temp);
    ratioM(1,iparti,iRatio) = mean(temp);
    ratioM(2,iparti,iRatio) = temp(190);
    ratioM(3,iparti,iRatio) = temp(10);
    fprintf('parti %d: BAS/AL=%f; Upper=%f; Lower=%f;\n',...
        iparti, ratioM(1,iparti,iRatio), ratioM(2,iparti,iRatio), ratioM(3,iparti,iRatio));

%     iRatio = 4; %cond AL/posterior-indep & order-dep
%     temp = parMC(3,iparti,:)./parMC(7,iparti,:);
%     temp = sort(temp);
%     ratioM(1,iparti,iRatio) = mean(temp);
%     ratioM(2,iparti,iRatio) = temp(190);
%     ratioM(3,iparti,iRatio) = temp(10);
%     fprintf('parti %d: AL/order-dep=%f; Upper=%f; Lower=%f;\n',...
%         iparti, ratioM(1,iparti,iRatio), ratioM(2,iparti,iRatio), ratioM(3,iparti,iRatio));
%     
%     iRatio = 5; %cond AL/posterior-dep & orer-indep
%     temp = parMC(3,iparti,:)./parMC(8,iparti,:);
%     temp = sort(temp);
%     ratioM(1,iparti,iRatio) = mean(temp);
%     ratioM(2,iparti,iRatio) = temp(190);
%     ratioM(3,iparti,iRatio) = temp(10);
%     fprintf('parti %d: AL/posterior-dep=%f; Upper=%f; Lower=%f;\n',...
%         iparti, ratioM(1,iparti,iRatio), ratioM(2,iparti,iRatio), ratioM(3,iparti,iRatio));
% 
%     iRatio = 5; %cond AL/posterior dep & order-dep
%     temp = parMC(3,iparti,:)./parMC(9,iparti,:);
%     temp = sort(temp);
%     ratioM(1,iparti,iRatio) = mean(temp);
%     ratioM(2,iparti,iRatio) = temp(190);
%     ratioM(3,iparti,iRatio) = temp(10);
%     fprintf('parti %d: AL/both-dep=%f; Upper=%f; Lower=%f;\n',...
%         iparti, ratioM(1,iparti,iRatio), ratioM(2,iparti,iRatio), ratioM(3,iparti,iRatio));
    
    fprintf('\n');
end
    
end
%%
% [SC1,~,~] = DATAFILE_Read('data\SC20130902ALS1.DAT');
% [SC2,~,~] = DATAFILE_Read('data\SC20130902ALS2.DAT');
% [SC3,~,~] = DATAFILE_Read('data\SC20130902ALS3.DAT');
% [SC4,~,~] = DATAFILE_Read('data\SC20130902ALS4.DAT');
% [SC5,~,~] = DATAFILE_Read('data\SC20130902ALS5.DAT');
% 
% [SC6,~,~] = DATAFILE_Read('data\SC20130903ALM2.DAT');
% [SC7,~,~] = DATAFILE_Read('data\SC20130903ALM3.DAT');
% [SC8,~,~] = DATAFILE_Read('data\SC20130903ALM4.DAT');
% [SC9,~,~] = DATAFILE_Read('data\SC20130904ALM5.DAT');
% [SC10,~,~] = DATAFILE_Read('data\SC20130904ALM6.DAT');
% [SC11,~,~] = DATAFILE_Read('data\SC20130904ALM8.DAT');
%  
% [EP1,~,~] = DATAFILE_Read('data\EP20130902ALS1.DAT');
% [EP2,~,~] = DATAFILE_Read('data\EP20130902ALS2.DAT');
% [EP3,~,~] = DATAFILE_Read('data\EP20130902ALS3.DAT');
% [EP4,~,~] = DATAFILE_Read('data\EP20130902ALS4.DAT');
% [EP5,~,~] = DATAFILE_Read('data\EP20130902ALS5.DAT');
% 
% [EP6,~,~] = DATAFILE_Read('data\EP20130903ALM1.DAT');
% [EP7,~,~] = DATAFILE_Read('data\EP20130903ALM2.DAT');
% [EP8,~,~] = DATAFILE_Read('data\EP20130903ALM3.DAT');
% [EP9,~,~] = DATAFILE_Read('data\EP20130904ALM4.DAT');
% [EP10,~,~] = DATAFILE_Read('data\EP20130904ALM6.DAT');
% [EP11,~,~] = DATAFILE_Read('data\EP20130904ALM7.DAT');
% 
% [BZ1,~,~] = DATAFILE_Read('data\BZ20130905ALS1.DAT');
% [BZ2,~,~] = DATAFILE_Read('data\BZ20130905ALS2.DAT');
% [BZ3,~,~] = DATAFILE_Read('data\BZ20130905ALS3.DAT');
% [BZ4,~,~] = DATAFILE_Read('data\BZ20130906ALS4.DAT');
% [BZ5,~,~] = DATAFILE_Read('data\BZ20130906ALS5.DAT');
% 
% [BZ6,~,~] = DATAFILE_Read('data\BZ20130906ALM1.DAT');
% [BZ7,~,~] = DATAFILE_Read('data\BZ20130906ALM2.DAT');
% [BZ8,~,~] = DATAFILE_Read('data\BZ20130906ALM3.DAT');
% [BZ9,~,~] = DATAFILE_Read('data\BZ20130907ALM5.DAT'); BZAL5.FrameData=[];
% [BZ10,~,~] = DATAFILE_Read('data\BZ20130907ALM6.DAT'); BZAL6.FrameData=[];
% [BZ11,~,~] = DATAFILE_Read('data\BZ20130907ALM7.DAT'); BZAL7.FrameData=[];
% 
% % DActive has 2 layers: parti x3; seesions
% for parti = 1:3
%     if (parti==1)       name = 'SC';
%     elseif (parti==2)   name = 'EP';
%     elseif (parti==3)   name = 'BZ';
%     end
%     for ss = 1:11
%         DActive(parti,ss) = eval([name, num2str(ss)]);
%     end
% end

% for parti = 4
%     for ss = 1:11
%         DActive(parti,ss) = a_Combine_Data(DActive(1,ss),DActive(2,ss),DActive(3,ss));
%     end
% end
%% test different functions

% for a = 1:0.1:5
%     xx = 1:0.1:25;
%     yy = 2./(1 + exp(-(xx-1)/exp(a)))-0.5;
%     plot(xx,yy); hold on;
% end




