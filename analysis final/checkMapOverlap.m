%function [paIm, shIm, svIm] = checkMapOverlap()
load('SCdata.mat', 'SCAL');

SCAL.RevealPosX(SCAL.RevealPosX==0) = nan;
SCAL.RevealPosY(SCAL.RevealPosY==0) = nan;

xDx = 1+(770-1)*(SCAL.RevealPosX/20+0.5);
xDy = 1-(770-1)*(SCAL.RevealPosY/20-0.5);
xDx = round(xDx);
xDy = round(xDy);

pa = SCAL.ImageID(:,2)==20;
sh = SCAL.ImageID(:,2)==6;
sv = SCAL.ImageID(:,2)==30;

paSubRow = xDx(pa,:);
paSubCol = xDy(pa,:);

shSubCol = xDy(sh,:);
shSubRow = xDx(sh,:);

svSubRow = xDx(sv,:);
svSubCol = xDy(sv,:);
 
paLinInd = sub2ind([770,770], paSubRow(:), paSubCol(:));
shLinInd = sub2ind([770,770], shSubRow(:), shSubCol(:));
svLinInd = sub2ind([770,770], svSubRow(:), svSubCol(:));

paLinInd(isnan(paLinInd))=[];
shLinInd(isnan(shLinInd))=[];
svLinInd(isnan(svLinInd))=[];

paIm = zeros(770);
for el = 1:numel(paLinInd)
    paIm(paLinInd(el)) = paIm(paLinInd(el))+1;
end
shIm = zeros(770);
for el = 1:numel(shLinInd)
    shIm(shLinInd(el)) = shIm(shLinInd(el))+1;
end
svIm = zeros(770);
for el = 1:numel(svLinInd)
    svIm(svLinInd(el)) = svIm(svLinInd(el))+1;
end

allIm = paIm + shIm + svIm;
allIm = allIm(:);
allIm(allIm==0)=[];
allIm = sort(allIm);
fracTwoOrMoreFixations = sum(allIm>2)/numel(allIm);

%end

