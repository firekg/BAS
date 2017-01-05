function [SIM] = simu_BASs0(sigp)
% for a quick a messy check of noisless BAS
%% SIM has 2 layers: participants x3 (SC,EP,BZ); conditions x3 (BAS, mBAS, pnBAS)
    biasw=zeros(26,1);
    %pv= [ lsigo,  phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, w];
    pv =[log(1),   nan,     nan,   nan,        nan,        nan,     0,  0.5, 1];
    %pw=[lbm,  km,  nm, bias];
    pw=[   0, nan, nan,    0]; % none of the last 3 are used here
    % make probe points z
    npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
    [x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
    dxx = xx(2)-xx(1);

    for i=1:1
        if(i==1)      
            load('SCdata.mat');  
            AL=SCAL;  
            %PLR=SCPLRbt;  
            DIMAL=a_DIM(SCAL);  
            %DIMPLR=a_DIM(SCPLRbt);
            pv(1) = log(sigp);  % from max likelihood fit
            disp(sprintf('Done loading.'));
        elseif(i==2)  
            load('EPdata.mat');  
            AL=EPAL;  
            PLR=EPPLRbt;  
            DIMAL=a_DIM(EPAL);  
            DIMPLR=a_DIM(EPPLRbt);
            pv(1) = log(1.2/10);  % from max likelihood fit
        elseif(i==3)  
            load('BZdata.mat');  
            AL=BZAL;  
            PLR=BZPLRbt;  
            DIMAL=a_DIM(BZAL);  
            DIMPLR=a_DIM(BZPLRbt);
            pv(1) = log(0.9/10); % from max likelihood fit
        end
        pw(1)=inf;     SIM(i,1) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,3);  %BAS
        %pw(1)=inf;     SIM(i,2) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %mBAS = BAS + motor noise
        %pw(1)=log(1);  SIM(i,3) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %nBAS = BAS + motor noise + selection noise
        disp(sprintf('parti %d SIM done.',i));
    end

end

%% These are the codes used: results in revMap0_017
% [SIM017] = simu_BASs0(0.17);
% [SIM0] = simu_BASs0(0);
% %%
% for i=1:3
%     for j=1:3
%         SIMd(i,j)=SIM017;
%     end
% end

