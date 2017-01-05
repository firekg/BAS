function [INFO] = info_BASs(SIM,deciSIM)
%% This is mainly to justify our usage of sig_p = 0.3 rather than 1.
% Want to show that the information gain is similar between these two types of revealings for various noise levels.
% They are mnn... OK similar?
% The issue is the revealings using sig_p=1 isn't robust at all!
% So, many other revealings would beat it....
% PROP PROP has 2 layers: noise leve x3 (1,0.3,0.1); conditions x3 (PLBprop, BASprop, deciBASprop)
% INFO has 4 layers: revealing number (1:25); parti x3; cond x3 (PLB, BAS, deciBAS); moments x2 (mean,SD)

[PROP] = prop_BASs(SIM,deciSIM);

INFO = nan(25,4,5,2);
for parti = 1:3
    for cond = 1:3
        m=PROP(parti, cond).mla;  
        M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  
        M= -M/log(0.5);
        INFO(:, parti, cond, 1) = nanmean(M);  
        INFO(:, parti, cond, 2) = sqrt(nanvar(M)./sum(~isnan(M)));        
    end
   %disp(sprintf('parti=%d;',parti));
end

end

%% PROP has 2 layers: noise leve x3 (1,0.3,0.1); conditions x3 (PLBprop, BASprop, deciBASprop)
function [PROP] = prop_BASs(SIM,deciSIM)
    Gmap=1;
    %pv= [ lsigo,  phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, w];
    pv =[ log(1),   nan,     nan,   nan,        nan,        nan,     0,  0.5, 1];
    %pw=[lbm,  km,  nm, bias];
    pw=[   0, nan, nan,    0]; % none of the last 3 are used here
    % make probe points z
    npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
    [x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
    dxx = xx(2)-xx(1);
    
    load('SCdata.mat');  
    PLB=SCPLBbt;  
    DIMPLB=a_DIM(PLB);  
    DIMAL=a_DIM(SCAL);  
    
    for i=1:3
        if(i==1)      
            pv(1) = log(1);
        elseif(i==2)  
            pv(1) = log(0.3);
        elseif(i==3)  
            pv(1) = log(0.1);
        end
        PROP(i,1) = a_Get_BALDscoreProp(PLB,DIMPLB,pv,pw,randn(1000,25),Gmap,z,2);  %PLB prop
        disp(sprintf('Done PLB.'));
        PROP(i,2) = a_Get_BALDscoreProp(SIM(1,1),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS for SCAL
        disp(sprintf('Done BAS.'));
        PROP(i,3) = a_Get_BALDscoreProp(deciSIM(1,1),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %deciBAS for SCAL
        disp(sprintf('Done deciBAS.'));
        disp(sprintf('%d PROP done.',i));
    end

end