function [SIM] = simu_MaxEnt()
%% SIM has 2 layers: participants x3 (SC,EP,BZ); conditions x1 (Max Ent)
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
            DIMAL=a_DIM(SCAL);  
            pv(1) = log(1.1);  % from max likelihood fit
        elseif(i==2)  
            load('EPdata.mat');  
            AL=EPAL;  
            DIMAL=a_DIM(EPAL);  
            pv(1) = log(1.2);  % from max likelihood fit
        elseif(i==3)  
            load('BZdata.mat');  
            AL=BZAL;  
            DIMAL=a_DIM(BZAL);  
            pv(1) = log(0.9); % from max likelihood fit
        end        
        pw(1)=inf;     SIM(i,1) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,5,[],biasw);  %Max Ent
        disp(sprintf('parti %d SIM done.',i));
    end

end