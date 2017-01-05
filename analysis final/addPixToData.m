%% copied from various places
function Dout = addPixToData(D)

DIM = a_DIM(D);
npix = 770;
ImageSize = 20;
revMax = 25;
yD = nan(D.Trials, revMax);

for trial = 1:D.Trials
    
    nr = D.MaxRevealingTrial(trial);
    xDx = 1+(npix-1)*(D.RevealPosX(trial,1:nr)/ImageSize+0.5);  % cm to pixel
    xDy = 1-(npix-1)*(D.RevealPosY(trial,1:nr)/ImageSize-0.5);
    xDx = round(xDx);
    xDy = round(xDy);
    xDx(xDx < 1) = 1;
    xDy(xDy < 1) = 1;
    xDx(xDx > npix) = npix;
    xDy(xDy > npix) = npix;
    
    yimage = DIM(trial).yimage;
    
    for rev = 1:nr
        yD(trial, rev) = yimage(xDy(rev), xDx(rev));  % x is column of image, y is row of image;
    end
    
end

Dout = D;
% to conform to the rest of SIMwPix
Dout.ml = [];
Dout.winner = [];
Dout.RevealZ = yD;

end
%%
% load('SCdata.mat', 'SCAL', 'SCPLRbt', 'SCPLBbt', 'SCPLaBbt');  
% load('EPdata.mat', 'EPAL', 'EPPLRbt', 'EPPLBbt', 'EPPLaBbt');  
% load('BZdata.mat', 'BZAL', 'BZPLRbt', 'BZPLBbt', 'BZPLaBbt');  

% SCALwPix = addPixToData(SCAL);
% SCPLRwPix = addPixToData(SCPLRbt);
% SCPLBwPix = addPixToData(SCPLBbt);
% SCPLaBwPix = addPixToData(SCPLaBbt);
% 
% EPALwPix = addPixToData(EPAL);
% EPPLRwPix = addPixToData(EPPLRbt);
% EPPLBwPix = addPixToData(EPPLBbt);
% EPPLaBwPix = addPixToData(EPPLaBbt);
% 
% BZALwPix = addPixToData(BZAL);
% BZPLRwPix = addPixToData(BZPLRbt);
% BZPLBwPix = addPixToData(BZPLBbt);
% BZPLaBwPix = addPixToData(BZPLaBbt);

% all the wPix data are saved as expDataWPix.mat
