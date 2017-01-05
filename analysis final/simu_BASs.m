function [SIM] = simu_BASs()
%% new version 2015-05-08
% SIM has 2 layers: participants x3 (SC,EP,BZ); conditions x3 (BAS, mBAS, pnBAS)

    biasw = zeros(26,1);    
    %pv= [  lsigo,    phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, window, prxpa, prCommon, shift];
    pv =[log(1.0),   nan,     nan,   nan,       log(0),        nan,     0,  1/2,      1,   1/3,      1/3,    12];
    %pw=[lbm,  km,  nm, bias];
    pw=[   0, nan, nan,    0]; % none of the last 3 are used here
    
    % make probe points z
    npix=770;
    ImageSize=20;
    n=110;
    xx = linspace(1,npix,n);
    yy = linspace(1,npix,n);
    [x1,x2]=meshgrid(xx,yy);
    z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
    dxx = xx(2)-xx(1);

    for i=1:3
        if(i==1)      
            load('SCdata.mat', 'SCAL');  
            AL=SCAL;  
            DIMAL=a_DIM(SCAL);  
            pv(1) = log(0.5);
            pv(6) = 0.34901995;
            pv(7) = 0.04471412;
            pv(12) = 16;
        elseif(i==2)  
            load('EPdata.mat', 'EPAL');  
            AL=EPAL;  
            DIMAL=a_DIM(EPAL);  
            pv(1) = log(0.5);
            pv(6) = 0.637376715;
            pv(7) = 0.120818875;
            pv(12) = 17;
        elseif(i==3)  
            load('BZdata.mat', 'BZAL');  
            AL=BZAL;  
            DIMAL=a_DIM(BZAL);  
            pv(1) = log(0.3);
            pv(6) = 0.374908933;
            pv(7) = 0.101003052;
            pv(12) = 15;
        end        
        for task = 1%:10
            
%             pw(1)=inf;
%             SIM(i,2, task) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %mBAS = BAS + motor noise
%             pw(1)=log(1);
%             SIM(i,3, task) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %nBAS = BAS + motor noise + selection noise
%             SIM(i,4, task) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,7);    % hardmax poster-dep & order-indep fixations
%             SIM(i,5, task) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,7.5);  % hardmax poster-dep & order-dep fixations            
% saved as analysisData-2015-06-16

%             pw(1) = inf;
%             pv(1) = log(0.17);
%             SIM(i,1) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw); %ideal BAS with prior bias
%             pv(12) = 0;
%             SIM(i,2) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %ideal BAS without prior bias
% saved in analysisPROP-priorBias-2015-07-02

% answering reviewer question-2015-09-04
            pw(1) = inf;
            SIM(i,2, task) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,5,[],biasw);  %mMaxEnt = MaxEnt + motor noise

            fprintf('parti %d, task %d SIM done.\n', i, task);
        end
    end

end

%% for comparing ideal BAS with and without prior bias - 2015-07-02
% pw(1) = inf;
% pv(1) = log(0.17);
% SIM(i,1) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw); %ideal BAS with prior bias
% pv(12) = 0;
% SIM(i,2) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %ideal BAS without prior bias
% Also need to manually change 2 things in other pieces of code
% 1. do not add perception noise to pixel values (a_Model Line 64)
% 2. do not add motor noise to selected location (a_Get_BALDscore Line 290)
% Remember to change them back!!!

%% old version
% function [SIM] = simu_BASs()
% %% SIM has 2 layers: participants x3 (SC,EP,BZ); conditions x3 (BAS, mBAS, pnBAS)
%     biasw=zeros(26,1);
%     %pv= [ lsigo,  phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, w];
%     pv =[log(1),   nan,     nan,   nan,        nan,        nan,     0,  0.5, 1];
%     %pw=[lbm,  km,  nm, bias];
%     pw=[   0, nan, nan,    0]; % none of the last 3 are used here
%     % make probe points z
%     npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
%     [x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
%     dxx = xx(2)-xx(1);
% 
%     for i=1:3
%         if(i==1)      
%             load('SCdata.mat', 'SCAL');  
%             AL=SCAL;  
%             DIMAL=a_DIM(SCAL);  
%             %pv(1) = log(1.1);  % from max likelihood fit
%             pv(1) = log(0.0687);  % from color difference exp
%         elseif(i==2)  
%             load('EPdata.mat', 'EPAL');  
%             AL=EPAL;  
%             DIMAL=a_DIM(EPAL);  
%             %pv(1) = log(1.2);  % from max likelihood fit
%             pv(1) = log(0.0687);  % from color difference exp
%         elseif(i==3)  
%             load('BZdata.mat', 'BZAL');  
%             AL=BZAL;  
%             DIMAL=a_DIM(BZAL);  
%             %pv(1) = log(0.9); % from max likelihood fit
%             pv(1) = log(0.0687);  % from color difference exp
%         end        
%         %pw(1)=inf;     SIM(i,1) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,3);  %BAS
%          pw(1)=inf;     SIM(i,2) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %mBAS = BAS + motor noise
%          pw(1)=log(1);  SIM(i,3) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %nBAS = BAS + motor noise + selection noise
%         %pw(1)=inf;     SIM(i,4) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,7);    % hardmax w-dep & order-indep fixations
%         %pw(1)=inf;     SIM(i,5) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,7.5);  % hardmax w-dep & order-dep fixations
% %         pw(1)=inf;     SIM(i,5) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,8);  % order-indep revMap prob match posterior
% %         pw(1)=inf;     SIM(i,6) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,9);  % order-indep saccadeMap hardmax posterior
% %         pw(1)=inf;     SIM(i,7) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,10); % order-indep saccadeMap prob match posterior
% %         pw(1)=inf;     SIM(i,8) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,11); % gungho hardmax posterior
% %         pw(1)=inf;     SIM(i,9) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,12); % gungho prob match posterior
%         disp(sprintf('parti %d SIM done.',i));
%     end
% 
% end