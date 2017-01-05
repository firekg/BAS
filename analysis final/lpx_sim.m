function [lpx] = lpx_sim()
%% px has 2 layers: participants x3 (SC,EP,BZ); conditions x3 (PLR, BAS, antiBAS)
    load('logNorRevMap.mat', 'logNorMap');
    for parti=1:3
        if(parti==1)      
            load('SCdata.mat', 'SCPLRbt', 'SCPLBbt', 'SCPLaBbt');  
            PLR=SCPLRbt;  
            PLB=SCPLBbt;  
            PLaB=SCPLaBbt;  
        elseif(parti==2)  
            load('EPdata.mat', 'EPPLRbt', 'EPPLBbt', 'EPPLaBbt');  
            PLR=EPPLRbt;  
            PLB=EPPLBbt;  
            PLaB=EPPLaBbt;  
        elseif(parti==3)  
            load('BZdata.mat', 'BZPLRbt', 'BZPLBbt', 'BZPLaBbt');  
            PLR=BZPLRbt;  
            PLB=BZPLBbt;  
            PLaB=BZPLaBbt;  
        end
        fprintf('Done loading parti %d.\n', parti);
        lpx(parti,1) = calLpx(PLR, logNorMap);  %random
        lpx(parti,2) = calLpx(PLB, logNorMap);  %BAS
        lpx(parti,3) = calLpx(PLaB, logNorMap);  %anti-BAS
        fprintf('Done px parti %d.\n', parti);
    end

end
%% copy stuff from a_Model
function lpx = calLpx(D, map)

npix = 770;
ImageSize = 20;
lpx.lik = nan(D.Trials, 25, 3);
lpx.mla = nan(D.Trials, 25);
lpx.Trials = D.Trials;

for trial = 1:D.Trials
    
    revMax = D.MaxRevealingTrial(trial);
    xDx = 1+(npix-1)*(D.RevealPosX(trial,1:revMax)/ImageSize+0.5);  % cm to pixel
    xDy = 1-(npix-1)*(D.RevealPosY(trial,1:revMax)/ImageSize-0.5);
    xDx = round(xDx);
    xDy = round(xDy);
    xDx(xDx < 1) = 1;
    xDy(xDy < 1) = 1;
    xDx(xDx > npix) = npix;
    xDy(xDy > npix) = npix;
    
    loglikPx = nan(revMax,3);
    for nr = 1:revMax
        for rev = 1:nr
            loglikPx(rev, 1) = map(xDy(rev), xDx(rev), 1, nr); % PA; x is column of image, y is row of image
            loglikPx(rev, 2) = map(xDy(rev), xDx(rev), 2, nr); % SH
            loglikPx(rev, 3) = map(xDy(rev), xDx(rev), 3, nr); % SV
        end
        logLikPx = nansum(loglikPx,1);
        lpx.lik(trial, nr, :) = logLikPx;
        lNor = log(sum( 0.5*exp(logLikPx(1)) + 0.25*exp(logLikPx(2)) + 0.25*exp(logLikPx(3)) )); %hardwiring unbiased prpa
        lpx.mla(trial, nr) = exp( log(0.5)+logLikPx(1)-lNor );
    end
    %fprintf('Done trial %d.\n', trial);
        
end

end

