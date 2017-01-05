function [PROP] = prop_progress(DALwPix)
%% PROP has 2 layers: participants x3 (SC,EP,BZ); conditions x5 (ALprop, PLRprop, BASprop, mBASprop, nBASprop)    
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
        
    for parti = 1:3 %3
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
        for session = 1:11
        D = DALwPix(parti, session);
            for iboot = 1:20
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
                for samp = 1:20       
                    PROP(parti, session, iboot, samp) = a_Get_BALDscoreProp(AFboot,[],pv,pw,randn(1000,25),Gmap,z,2);
                end
            end    
        fprintf('Done parti %d, session %d.\n', parti, session);
        end
    end

end
