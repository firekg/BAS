%%

function [posX, posY] = transformPosInd2Cm(posRow, posCol)
    expConst = defineExperimentConstants();
    [pixelsPerDim, ~, imageLengthCm] = unpackStruct(expConst);
    posX = (posCol-1)/(pixelsPerDim-1)*imageLengthCm - imageLengthCm/2;
    posY = (1-posRow)/(pixelsPerDim-1)*imageLengthCm + imageLengthCm/2;
end